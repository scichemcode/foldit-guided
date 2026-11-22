import recipe from './recipe.js';
import demo from './demo.js';
import * as store from './store.js';
import { run } from './engine.js';
import { findCaptures } from './physics/capture.js';
import { parsePdb } from './physics/pdb.js';
import { exportSteps, importSteps } from './transfer.js';
import { canFollow, follow, pickLog, rememberedLog } from './watch.js';
import { explainLog, markDone, setUpGuide, updateGuide } from './guide.js';
import { NOISE, summarise } from './summary.js';

const $ = selector => document.querySelector(selector);
const fmt = (v, d = 2) => (v == null ? '–' : v.toFixed(d));
const signed = (v, d = 2) => (v == null ? '–' : (v > 0 ? '+' : '') + v.toFixed(d));

// Names, labels and file contents are untrusted: they only ever become text nodes.
function el(tag, props = {}, ...kids) {
  const node = Object.assign(document.createElement(tag), props);
  node.append(...kids.filter(k => k != null && k !== false));
  return node;
}
function svg(tag, attrs, text) {
  const node = document.createElementNS('http://www.w3.org/2000/svg', tag);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
  if (text != null) node.textContent = text;
  return node;
}

let steps = [], terms = [], challengeId = null, stepKey = null, baselineKey = null, latest = true, drawing = 0;
let guideOpen = null; // null until we know whether this is a first visit
let lastOrder = 0;
const nextOrder = () => (lastOrder = Math.max(lastOrder + 1, Date.now()));
const say = text => { $('#message').textContent = text; };
const inChallenge = () => steps.filter(s => s.challenge === challengeId);

function challenges() {
  const map = new Map();
  for (const s of steps) {
    const c = map.get(s.challenge) ?? { id: s.challenge, count: 0, best: null };
    c.name = s.name || c.name;
    c.count++;
    if (s.score != null) c.best = Math.max(c.best ?? -Infinity, s.score);
    c.order = s.order;
    map.set(s.challenge, c);
  }
  return [...map.values()].sort((a, b) => b.order - a.order);
}

// The previous analysed step on the same track of the same challenge.
function previousSteps(list) {
  const last = new Map(), previous = new Map();
  for (const s of list) {
    if (s.status !== 'done') continue;
    if (last.has(s.track)) previous.set(s.key, last.get(s.track));
    last.set(s.track, s);
  }
  return previous;
}

const shapeOf = s => s.residues.map(r => `${r.name}${r.id}/${r.chain}`).join() + s.x.length;
const energyChange = (from, to) => (from && shapeOf(from) === shapeOf(to) ? to.energy.total - from.energy.total : null);
const scoreChange = (from, to) => (from?.score != null && to.score != null ? to.score - from.score : null);

function verdict(ds, de) {
  if (ds == null || de == null || Math.abs(ds) < 0.005 || Math.abs(de) < NOISE) return '–';
  return (ds > 0) === (de < 0) ? 'agree' : 'disagree';
}

async function reload() {
  steps = await store.allSteps();
  lastOrder = Math.max(lastOrder, ...steps.map(s => s.order));
  const all = challenges();
  if (!all.some(c => c.id === challengeId)) {
    challengeId = all[0]?.id ?? null;
    latest = true;
  }
  await render();
}

async function addAndAnalyse(list, options) {
  const added = await store.addSteps(list, options);
  if (!added) return 0;
  challengeId = list.at(-1).challenge;
  guideOpen = false;
  baselineKey = null;
  latest = true;
  await reload();
  processQueue();
  return added;
}

async function addText(text) {
  const list = findCaptures(text).map(c => ({
    key: c.key, challenge: `puzzle:${c.puzzle}`, name: c.name || c.puzzle, description: c.description, track: c.track,
    move: c.move, score: c.score, parts: c.parts, notes: c.warnings, raw: c.raw, status: 'pending', order: nextOrder(),
  }));
  return list.length ? addAndAnalyse(list) : 0;
}

let working = false, queuedAgain = false;
async function processQueue() {
  if (working) { queuedAgain = true; return; }
  working = true;
  const drain = async () => {
    for (;;) {
      const step = (await store.allSteps()).find(s => s.status === 'pending');
      if (!step) break;
      let done;
      try {
        const out = await run('analyse', { raw: step.raw, residues: step.residues }, iteration => {
          $('#progress').textContent = `Analysing ${step.move || 'capture'}: relaxing hydrogens, iteration ${iteration}`;
        });
        done = { ...step, raw: undefined, status: 'done', residues: out.residues, x: out.x, energy: out.energy, gaps: out.warnings };
      } catch (error) {
        done = { ...step, status: 'failed', error: error.message };
      }
      await store.updateStep(done);
      await reload();
    }
  };
  try {
    // The lock stops two open tabs from analysing the same step.
    do {
      queuedAgain = false;
      if (navigator.locks) await navigator.locks.request('foldit-energy-queue', drain);
      else await drain();
    } while (queuedAgain);
  } catch (error) {
    say(`Analysis stopped: ${error.message}`);
  } finally {
    working = false;
    $('#progress').textContent = '';
  }
}

function drawChallenges() {
  $('#challenges').replaceChildren(...challenges().map(c => {
    const button = el('button', { className: 'challenge', onclick: () => {
      challengeId = c.id; baselineKey = null; latest = true; guideOpen = false; say(''); render();
    } }, c.name || c.id, el('small', {},
      `${c.count} step${c.count === 1 ? '' : 's'}${c.best == null ? '' : ` · best score ${fmt(c.best, 1)}`}`));
    button.setAttribute('aria-current', !guideOpen && c.id === challengeId);
    return button;
  }));
  $('#export').disabled = !steps.length;
  $('#no-challenges').hidden = steps.length > 0;
}

function selectStep(key) {
  stepKey = key;
  baselineKey = null;
  latest = key === inChallenge().at(-1).key;
  render().then(() => $('#step-title')?.scrollIntoView({ block: 'nearest' }));
}

function history(title, list, value) {
  const points = list.map((s, i) => ({ s, i, v: value(s) })).filter(p => p.v != null);
  if (!points.length) return null;
  const lo = Math.min(...points.map(p => p.v)), hi = Math.max(...points.map(p => p.v)), span = hi - lo || 1;
  const x = i => (list.length === 1 ? 190 : 64 + (i * 284) / (list.length - 1));
  const y = v => 122 - ((v - lo) / span) * 108;
  const plot = svg('svg', { viewBox: '0 0 360 148', role: 'img', 'aria-label': title });
  for (const v of new Set([lo, hi])) {
    plot.append(svg('line', { class: 'grid', x1: 64, x2: 348, y1: y(v), y2: y(v) }),
      svg('text', { x: 56, y: y(v) + 4, 'text-anchor': 'end' }, fmt(v, 1)));
  }
  for (const track of new Set(points.map(p => p.s.track))) { // tracks are separate histories
    plot.append(svg('polyline', { class: 'line',
      points: points.filter(p => p.s.track === track).map(p => `${x(p.i)},${y(p.v)}`).join(' ') }));
  }
  const readout = el('div', { className: 'readout' });
  for (const { s, i, v } of points) {
    plot.append(svg('circle', { class: s.key === stepKey ? 'dot selected' : 'dot', cx: x(i), cy: y(v), r: 5 }));
    const hit = svg('circle', { class: 'hit', cx: x(i), cy: y(v), r: 14 });
    hit.onmouseenter = () => { readout.textContent = `Step ${i + 1} ${s.move || ''}: ${fmt(v)}`; };
    hit.onmouseleave = () => { readout.textContent = ''; };
    hit.onclick = () => selectStep(s.key);
    plot.append(hit);
  }
  plot.append(svg('text', { x: 206, y: 146, 'text-anchor': 'middle' }, `steps 1–${list.length}`));
  return el('figure', {}, el('figcaption', {}, title), plot, readout);
}

function termTable(delta) {
  const max = Math.max(1e-9, ...terms.map(t => Math.abs(delta.terms[t.key])));
  return el('div', { className: 'scroll' }, el('table', { id: 'terms' },
    el('thead', {}, el('tr', {}, el('th', {}, 'What changed'), el('th', {}, '← better | worse →'),
      el('th', {}, 'kcal/mol'), el('th', { className: 'l' }, 'Closest Foldit part'))),
    el('tbody', {}, ...terms.map(t => {
      const v = delta.terms[t.key], width = (Math.abs(v) / max) * 50;
      const fill = el('div', { className: v < 0 ? 'fill down' : 'fill up' });
      fill.style.width = `${width}%`;
      fill.style.marginLeft = `${v < 0 ? 50 - width : 50}%`;
      return el('tr', { title: t.equation }, el('td', {}, t.plain, el('small', {}, t.label)),
        el('td', { className: 'bar' }, el('div', { className: 'track' }, fill)),
        el('td', {}, signed(v)), el('td', { className: 'l' }, t.foldit));
    }))));
}

function explanation(delta) {
  if (!delta.comparable) {
    return el('p', { className: 'warn' }, delta.error ?? 'Atoms or bonds changed (mutation, disulfide, insertion or deletion). '
      + 'Energies of different molecules have different zeros, so no ΔE is given. The next capture will be comparable again.');
  }
  return el('section', { id: 'report' },
    delta.heading && el('h3', {}, delta.heading),
    delta.equation && el('p', { className: 'equation' }, delta.equation),
    ...delta.causes.map(c => el('div', { className: 'cause' },
      el('strong', {}, c.title), el('ul', {}, ...c.notes.map(n => el('li', {}, n))))),
    delta.hint && el('p', {}, delta.hint),
    delta.gain && el('p', {}, delta.gain));
}

async function render() {
  const turn = ++drawing;
  const list = inChallenge();
  guideOpen ??= !list.length;
  if (!list.length) guideOpen = true;
  drawChallenges();
  updateGuide(steps);
  $('#guide').hidden = !guideOpen;
  $('#view').hidden = guideOpen;
  $('#show-guide').setAttribute('aria-current', guideOpen);
  if (guideOpen) return;

  if (latest || !list.some(s => s.key === stepKey)) stepKey = list.at(-1).key;
  const step = list.find(s => s.key === stepKey);
  const previous = previousSteps(list);
  const base = step.status === 'done'
    ? list.find(s => s.key === baselineKey && s.status === 'done' && s !== step) ?? previous.get(step.key) : null;

  let delta = null;
  if (base) {
    try {
      const pick = ({ residues, x, energy }) => ({ residues, x, energy });
      delta = await run('compare', { before: pick(base), after: pick(step) });
    } catch (error) {
      delta = { comparable: false, error: error.message };
    }
    if (turn !== drawing) return; // a newer render has started
  }

  const changes = list.map(s => {
    const from = s.status === 'done' ? previous.get(s.key) : null;
    const ds = scoreChange(from, s), de = from ? energyChange(from, s) : null;
    return { ds, de, says: verdict(ds, de) };
  });
  const judged = changes.filter(c => c.says !== '–');
  const tile = (value, label) => el('div', {}, el('b', {}, value), el('span', {}, label));

  const picker = el('select', { id: 'baseline', onchange: e => { baselineKey = e.target.value; render(); } },
    ...list.filter(s => s !== step && s.status === 'done').map(s =>
      el('option', { value: s.key, selected: s === base }, `${list.indexOf(s) + 1}. ${s.move || 'capture'}`)));
  if (!base) picker.prepend(el('option', { selected: true, disabled: true }, 'choose a step'));

  const before = new Map(base?.parts ?? []);
  const parts = (base ? step.parts ?? [] : []).filter(([k]) => before.has(k)).map(([k, v]) => [k, v - before.get(k)])
    .filter(([, v]) => Math.abs(v) >= 0.005).sort((a, b) => Math.abs(b[1]) - Math.abs(a[1]));

  const remove = () => {
    if (confirm('Delete this challenge and its history from this browser?')) store.deleteSteps(list).then(reload);
  };
  const removeStep = () => {
    if (confirm('Delete this step from this browser?')) store.deleteSteps([step]).then(reload);
  };

  const ds = scoreChange(base, step);
  $('#view').replaceChildren(...[
    el('div', { className: 'title' }, el('h2', {}, list[0].name || challengeId),
      el('button', { onclick: remove }, 'Delete challenge')),
    nextAction(list),
    el('section', { className: 'result' },
      el('div', { className: 'title' },
        el('h3', { id: 'step-title' }, `Step ${list.indexOf(step) + 1}: ${step.move || 'capture'}`),
        el('button', { onclick: removeStep }, 'Delete step')),
      ...resultCard(step, base, delta, ds),
      ...[...(step.notes ?? []), ...(step.gaps ?? [])].map(note => el('p', { className: 'warn' }, note)),
      step.status === 'done' && list.length > 1 && el('label', { className: 'inline' }, 'Compared with ', picker)),

    el('h3', {}, 'Your captures'),
    el('p', { className: 'muted' }, 'Click a row to see why that move helped or hurt.'),
    el('div', { className: 'scroll' }, el('table', { id: 'steps' },
      el('thead', {}, el('tr', {}, ...['Step', 'Track', 'Foldit score', 'Score change', 'Energy', 'Energy change', 'Foldit and physics'].map(h => el('th', {}, h)))),
      el('tbody', {}, ...list.map((s, i) => el('tr',
        { className: s === step ? 'pick selected' : 'pick', onclick: () => selectStep(s.key) },
        el('td', {}, el('button', {}, `${i + 1}. ${s.move || 'capture'}`)), el('td', {}, s.track || '–'),
        el('td', {}, fmt(s.score)), el('td', {}, signed(changes[i].ds)),
        el('td', {}, s.status === 'done' ? fmt(s.energy.total) : s.status),
        el('td', {}, signed(changes[i].de)), el('td', {}, changes[i].says)))))),

    delta?.comparable && el('h3', {}, 'What changed in this move'),
    delta?.comparable && el('p', { className: 'muted' }, 'Red bars hurt and blue bars helped. The longer the bar, the more it mattered.'),
    delta?.comparable && termTable(delta),
    delta?.comparable && fold('physics', 'Show the atoms and equations', explanation(delta)),
    parts.length && el('section', {}, el('h3', {}, 'What Foldit’s own score parts did'),
      el('table', {}, el('tbody', {}, ...parts.map(([k, v]) => el('tr', {}, el('td', {}, k), el('td', {}, signed(v))))))),

    el('h3', {}, 'Your progress'),
    el('div', { className: 'tiles' }, tile(String(list.length), 'captures'),
      tile(judged.length ? `${judged.filter(c => c.says === 'agree').length} of ${judged.length}` : '–', 'moves where Foldit and physics agree')),
    el('div', { className: 'charts' }, ...[
      history('Foldit score (higher is better)', list, s => s.score),
      history('Amber energy, kcal/mol (lower is better)', list, s => s.energy?.total),
    ].filter(Boolean)),
    fold('glossary', 'What do these words mean?', glossary()),
    list[0].description && fold('description', 'About this puzzle', el('p', {}, list[0].description)),
  ].filter(Boolean));
}

const unfolded = new Set(); // the view is rebuilt often; remember what the reader opened
function fold(id, title, content) {
  return el('details', { id, open: unfolded.has(id),
    ontoggle: event => (event.target.open ? unfolded.add(id) : unfolded.delete(id)) },
  el('summary', {}, title), content);
}

function glossary() {
  return el('dl', {}, ...[
    ['Capture', 'One snapshot of your protein, taken by the recipe.'],
    ['Baseline', 'The capture that a move is compared with. Normally the one just before it.'],
    ['Energy', 'How strained the protein is according to standard physics (Amber). Lower is better. Only the changes mean something, not the number itself.'],
    ['kcal/mol', 'The unit of energy. A typical hydrogen bond is worth 1 to 5 kcal/mol.'],
    ['LEU23', 'An amino acid and its segment number: leucine at segment 23.'],
    ['Agree or disagree', 'Whether Foldit’s score and the physics both call the move good, or both bad.'],
  ].flatMap(([word, meaning]) => [el('dt', {}, word), el('dd', {}, meaning)]));
}

function nextAction(list) {
  const id = list[0].challenge, last = list.at(-1);
  if (id === 'demo') {
    return el('p', { className: 'next' }, 'This is an example. Step 2 pushed atoms into each other and step 3 repaired it: click them below. ',
      el('button', { onclick: openGuide }, 'Now set up Foldit'));
  }
  if (!id.startsWith('puzzle:') || last.status === 'pending') return null;
  const again = canFollow ? '' : ' Then come back and choose the script log again.';
  const text = last.status === 'failed' ? 'That capture failed. Read the reason below, then run the recipe again.'
    : list.length === 1 ? 'Starting point saved. Now make one move in Foldit and run the recipe again.'
    : 'Next: make a move in Foldit and run the recipe again.';
  return el('p', { className: 'next' }, el('span', {}, text + again + ' '),
    !canFollow && el('button', { onclick: () => $('#log-file').click() }, 'Choose the script log'));
}

function resultCard(step, base, delta, ds) {
  if (step.status === 'pending') {
    return [el('p', { className: 'headline' }, 'Working out the energy…'),
      el('p', { className: 'muted' }, 'From a second to a few minutes, depending on the protein’s size.')];
  }
  if (step.status === 'failed') {
    return [el('p', { className: 'headline' }, 'This capture could not be analysed'), el('p', { className: 'warn' }, step.error)];
  }
  if (!base) {
    return [el('p', { className: 'headline' }, 'This is a starting point'),
      el('p', {}, 'There is nothing earlier on this track to compare it with. The explanation appears with the next capture.')];
  }
  if (!delta.comparable) return [el('p', { className: 'headline' }, 'This move can’t be compared'), explanation(delta)];
  const { tone, headline, lines } = summarise(delta, ds, terms);
  return [el('p', { className: `headline ${tone}` }, headline),
    el('dl', { className: 'summary' }, ...lines.flatMap(([label, text]) => [el('dt', {}, label), el('dd', {}, text)]))];
}

function openGuide() {
  guideOpen = true;
  say('');
  render();
  scrollTo(0, 0);
}

function download(name, text, type) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  el('a', { href: url, download: name }).click();
  setTimeout(() => URL.revokeObjectURL(url), 5000);
}

function guarded(action) {
  return async event => {
    try {
      await action(event);
    } catch (error) {
      if (error.name !== 'AbortError') say(error.message); // AbortError: picker closed
    }
  };
}

async function readLog(text) {
  const added = await addText(text);
  say(added ? `Received ${added} new capture${added === 1 ? '' : 's'}.` : explainLog(text));
}

function startFollowing(handle) {
  follow(handle, readLog, (text, expired = false) => {
    $('#connection').textContent = text;
    $('#reconnect').hidden = !expired;
    $('#reconnect').onclick = guarded(async () => {
      if ((await handle.requestPermission({ mode: 'read' })) === 'granted') startFollowing(handle);
    });
  });
}

function bind() {
  $('#recipe').value = recipe;
  $('#copy').onclick = guarded(async () => {
    try {
      await navigator.clipboard.writeText(recipe);
      say('Recipe copied. Save it as a plain text file ending in .lua, then use Load local recipe in Foldit.');
    } catch {
      $('#recipe-box').open = true;
      $('#recipe').select();
      say('Copy the selected text.');
    }
  });
  $('#save-recipe').onclick = () => {
    download('foldit-energy-capture.lua', recipe, 'application/octet-stream');
    markDone('downloaded');
    say('Recipe downloaded. Next is step 2, in Foldit.');
  };
  $('#captured').onclick = () => {
    markDone('captured');
    say('Good. Next is step 3: show this page where the script log is.');
  };
  $('#show-guide').onclick = openGuide;
  for (const button of document.querySelectorAll('[data-copy]')) {
    button.onclick = guarded(async () => {
      try {
        await navigator.clipboard.writeText(button.dataset.copy);
        say('Location copied.');
      } catch {
        say('Copying is blocked here: select the location and copy it by hand.');
      }
    });
  }
  $('#follow').onclick = guarded(async () => startFollowing(await pickLog()));

  $('#log-file').onchange = guarded(async event => {
    let added = 0, problem = '';
    for (const file of event.target.files) {
      const text = await file.text(), count = await addText(text);
      added += count;
      if (!count) problem = explainLog(text);
    }
    event.target.value = ''; // so the same log can be chosen again
    say(added ? `Received ${added} new capture${added === 1 ? '' : 's'}.` : problem);
  });

  $('#paste').onsubmit = guarded(async event => {
    event.preventDefault();
    const text = event.target.elements.text.value, added = await addText(text);
    say(added ? `Received ${added} new capture${added === 1 ? '' : 's'}.` : explainLog(text));
    if (added) event.target.reset();
  });

  $('#upload').onsubmit = guarded(async event => {
    event.preventDefault();
    const { challenge, move, file } = event.target.elements;
    const name = challenge.value.trim();
    if (!name) throw new Error('Give the challenge a name.');
    const { residues, warnings } = parsePdb(await file.files[0].text());
    const order = nextOrder();
    await addAndAnalyse([{ key: `pdb:${order}:${Math.random()}`, challenge: `upload:${name}`, name, track: '',
      move: move.value.trim(), score: null, parts: [], notes: warnings, residues, status: 'pending', order }]);
    say('PDB added.');
  });

  const showExample = guarded(async () => {
    scrollTo(0, 0);
    say('');
    const existing = steps.filter(s => s.challenge === 'demo');
    if (existing.length) {
      challengeId = 'demo'; guideOpen = false; baselineKey = null; latest = false;
      stepKey = (existing[1] ?? existing[0]).key; // the interesting one
      return render();
    }
    const stamp = nextOrder();
    await addAndAnalyse(demo.map(({ move, residues }, i) => ({ key: `demo:${stamp}:${i}`, challenge: 'demo',
      name: 'Example: synthetic peptide', description: 'A test peptide, not a Foldit puzzle, so it has no Foldit scores.',
      track: '', move, score: null, parts: [], notes: [], residues, status: 'pending', order: nextOrder() })));
    latest = false;
    stepKey = `demo:${stamp}:1`;
    await render();
  });
  $('#demo').onclick = showExample;
  $('#demo-guide').onclick = showExample;

  $('#export').onclick = () => {
    if (confirm('This file contains your solutions. Do not share it while those puzzles are active. Save it?')) {
      download('foldit-energy-history.json', exportSteps(steps), 'application/json');
    }
  };
  $('#import').onchange = guarded(async event => {
    const file = event.target.files[0];
    const added = await addAndAnalyse(importSteps(await file.text()), { restore: true });
    event.target.value = '';
    say(`Imported ${added} step(s).`);
  });
}

async function boot() {
  bind();
  try {
    await setUpGuide();
    terms = await run('terms');
    await reload();
  } catch (error) {
    say(`Cannot start: ${error.message}. Private windows may block browser storage or workers.`);
    return;
  }
  processQueue();
  const saved = await rememberedLog();
  if (saved?.granted) startFollowing(saved.handle);
  else if (saved) {
    $('#connection').textContent = `Last used: ${saved.handle.name}`;
    $('#reconnect').hidden = false;
    $('#reconnect').onclick = guarded(async () => {
      if ((await saved.handle.requestPermission({ mode: 'read' })) === 'granted') startFollowing(saved.handle);
    });
  }
  document.body.dataset.ready = 'true';
}

boot();
