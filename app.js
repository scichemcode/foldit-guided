import recipe from './recipe.js';
import demo from './demo.js';
import * as store from './store.js';
import { run } from './engine.js';
import { findCaptures } from './physics/capture.js';
import { parsePdb } from './physics/pdb.js';
import { exportSteps, importSteps } from './transfer.js';
import { canFollow, follow, pickLog, rememberedLog } from './watch.js';

const $ = selector => document.querySelector(selector);
const fmt = (v, d = 2) => (v == null ? '–' : v.toFixed(d));
const signed = (v, d = 2) => (v == null ? '–' : (v > 0 ? '+' : '') + v.toFixed(d));
const NOISE = 0.05; // kcal/mol; smaller changes get no verdict

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
      challengeId = c.id; baselineKey = null; latest = true; render();
    } }, c.name || c.id, el('small', {},
      `${c.count} step${c.count === 1 ? '' : 's'}${c.best == null ? '' : ` · best score ${fmt(c.best, 1)}`}`));
    button.setAttribute('aria-current', c.id === challengeId);
    return button;
  }));
  $('#export').disabled = !steps.length;
}

function selectStep(key) {
  stepKey = key;
  baselineKey = null;
  latest = key === inChallenge().at(-1).key;
  render();
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
    el('thead', {}, el('tr', {}, el('th', {}, 'Term'), el('th', {}, '← downhill | uphill →'),
      el('th', {}, 'ΔE kcal/mol'), el('th', { className: 'l' }, 'Nearest Foldit part'))),
    el('tbody', {}, ...terms.map(t => {
      const v = delta.terms[t.key], width = (Math.abs(v) / max) * 50;
      const fill = el('div', { className: v < 0 ? 'fill down' : 'fill up' });
      fill.style.width = `${width}%`;
      fill.style.marginLeft = `${v < 0 ? 50 - width : 50}%`;
      return el('tr', { title: t.equation }, el('td', {}, t.label),
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
  drawChallenges();
  const list = inChallenge();
  $('#intro').hidden = list.length > 0;
  $('#view').hidden = !list.length;
  if (!list.length) return;

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

  $('#view').replaceChildren(...[
    el('div', { className: 'title' }, el('h2', {}, list[0].name || challengeId),
      el('button', { onclick: remove }, 'Delete challenge')),
    list[0].description && el('p', { className: 'muted' }, list[0].description),
    el('div', { className: 'tiles' }, tile(String(list.length), 'steps'),
      tile(judged.length ? `${judged.filter(c => c.says === 'agree').length} of ${judged.length}` : '–', 'moves where Foldit and Amber agree')),
    el('div', { className: 'charts' }, ...[
      history('Foldit score (higher is better)', list, s => s.score),
      history('Amber energy, kcal/mol (lower is better)', list, s => s.energy?.total),
    ].filter(Boolean)),
    el('div', { className: 'scroll' }, el('table', { id: 'steps' },
      el('thead', {}, el('tr', {}, ...['Step', 'Track', 'Score', 'ΔScore', 'E', 'ΔE', 'Verdict'].map(h => el('th', {}, h)))),
      el('tbody', {}, ...list.map((s, i) => el('tr',
        { className: s === step ? 'pick selected' : 'pick', onclick: () => selectStep(s.key) },
        el('td', {}, el('button', {}, `${i + 1}. ${s.move || 'capture'}`)), el('td', {}, s.track || '–'),
        el('td', {}, fmt(s.score)), el('td', {}, signed(changes[i].ds)),
        el('td', {}, s.status === 'done' ? fmt(s.energy.total) : s.status),
        el('td', {}, signed(changes[i].de)), el('td', {}, changes[i].says)))))),
    el('div', { className: 'title' }, el('h2', { id: 'step-title' }, `Step ${list.indexOf(step) + 1}: ${step.move || 'capture'}`),
      el('button', { onclick: removeStep }, 'Delete step')),
    step.error && el('p', { className: 'warn' }, step.error),
    ...[...(step.notes ?? []), ...(step.gaps ?? [])].map(note => el('p', { className: 'warn' }, note)),
    step.status === 'pending' && el('p', {}, 'Analysing…'),
    step.status === 'done' && list.length > 1 && el('label', {}, 'Compared with ', picker),
    step.status === 'done' && !base && el('p', { className: 'muted' }, 'Baseline: the next capture will be compared with it.'),
    delta?.comparable && el('div', { className: 'tiles' }, tile(signed(delta.total), 'ΔE, kcal/mol'),
      tile(signed(scoreChange(base, step)), 'ΔFoldit score'), tile(verdict(scoreChange(base, step), delta.total), 'verdict')),
    delta?.comparable && termTable(delta),
    delta && explanation(delta),
    parts.length && el('section', {}, el('h3', {}, 'What Foldit’s own score parts did'),
      el('table', {}, el('tbody', {}, ...parts.map(([k, v]) => el('tr', {}, el('td', {}, k), el('td', {}, signed(v))))))),
  ].filter(Boolean));
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

function startFollowing(handle) {
  follow(handle, addText, (text, expired = false) => {
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
      say('Recipe copied. Paste it into a new Lua recipe in Foldit.');
    } catch {
      $('#recipe-box').open = true;
      $('#recipe').select();
      say('Copy the selected text.');
    }
  });
  $('#save-recipe').onclick = () => download('capture.lua', recipe, 'text/plain');

  $('#follow').hidden = !canFollow;
  $('#no-follow').hidden = canFollow;
  $('#follow').onclick = guarded(async () => startFollowing(await pickLog()));

  $('#log-file').onchange = guarded(async event => {
    let added = 0;
    for (const file of event.target.files) added += await addText(await file.text());
    event.target.value = ''; // so the same log can be chosen again
    say(added ? `Added ${added} capture(s).` : 'No new complete capture in that file.');
  });

  $('#paste').onsubmit = guarded(async event => {
    event.preventDefault();
    const added = await addText(event.target.elements.text.value);
    say(added ? `Added ${added} capture(s).` : 'No new complete capture in that text.');
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

  $('#demo').onclick = guarded(async () => {
    const existing = steps.find(s => s.challenge === 'demo');
    if (existing) {
      challengeId = 'demo'; latest = true;
      return render();
    }
    const stamp = nextOrder();
    await addAndAnalyse(demo.map(({ move, residues }, i) => ({ key: `demo:${stamp}:${i}`, challenge: 'demo',
      name: 'Example: synthetic peptide', description: 'A test peptide, not a Foldit puzzle, so it has no Foldit scores.',
      track: '', move, score: null, parts: [], notes: [], residues, status: 'pending', order: nextOrder() })));
  });

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
