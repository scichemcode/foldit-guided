// The step-by-step guide: which steps are done, and which one to show next.
import { getSetting, setSetting } from './store.js';
import { canFollow } from './watch.js';

const $ = selector => document.querySelector(selector);
const progress = { downloaded: false, captured: false };
let current = null, records = [];

export function guessSystem() {
  const platform = navigator.userAgentData?.platform || navigator.platform || '';
  if (/mac/i.test(platform)) return 'mac';
  return /linux/i.test(platform) ? 'linux' : 'windows';
}

export async function setUpGuide() {
  for (const node of document.querySelectorAll('.can-follow')) node.hidden = !canFollow;
  for (const node of document.querySelectorAll('.cannot-follow')) node.hidden = canFollow;

  const system = $('#system');
  const showSystem = () => {
    for (const node of document.querySelectorAll('[data-os]')) node.hidden = node.dataset.os !== system.value;
  };
  system.value = guessSystem();
  system.onchange = showSystem;
  showSystem();

  Object.assign(progress, await getSetting('guide').catch(() => null));
}

export function markDone(name) {
  progress[name] = true;
  updateGuide(records);
  setSetting('guide', { ...progress }).catch(() => {}); // the guide works even if it cannot be remembered
}

// Progress follows the player's real captures; the example and uploads don't count.
export function updateGuide(steps) {
  records = steps;
  const perPuzzle = new Map();
  for (const step of steps) {
    if (!step.challenge.startsWith('puzzle:') || step.status === 'failed') continue;
    perPuzzle.set(step.challenge, (perPuzzle.get(step.challenge) ?? 0) + 1);
  }
  const received = perPuzzle.size > 0, compared = [...perPuzzle.values()].some(count => count >= 2);
  const done = [progress.downloaded || progress.captured || received, progress.captured || received, received, compared];
  const next = done.indexOf(false) + 1; // 0: all done

  for (const step of document.querySelectorAll('.step')) {
    const number = Number(step.dataset.step), finished = done[number - 1];
    step.classList.toggle('done', finished);
    step.classList.toggle('current', number === next);
    step.querySelector('.number').textContent = finished ? '✓' : String(number);
    step.querySelector('.state').textContent = finished ? 'Done' : number === next ? 'Do this now' : '';
    if (next !== current) step.open = number === next; // only when progress changes, so reading ahead works
  }
  current = next;
  $('#step-count').textContent = next ? `You are on step ${next} of 4.` : 'Setup is done. Click any step to read it again.';
}

// What to tell the player about a log that added no new capture.
export function explainLog(text) {
  if (!text.trim()) return 'That file is empty. Run the capture recipe in Foldit first.';
  const started = text.lastIndexOf('FEX1|BEGIN'), ended = text.lastIndexOf('FEX1|END');
  const stopped = [...text.matchAll(/capture stopped: ([^\r\n<]*)/g)].at(-1);
  if (stopped && stopped.index > ended) {
    return `The recipe stopped before finishing (${stopped[1].slice(0, 160)}). Run it again.`;
  }
  if (started > ended) {
    const percent = [...text.slice(started).matchAll(/Capturing: (\d+)%/g)].at(-1);
    return `A capture is under way${percent ? ` (${percent[1]}%)` : ''}. Wait for “Capture finished” in Foldit.`;
  }
  if (ended >= 0) return 'No new capture: this one was already added. Make a move, then run the recipe again.';
  return 'No capture in that file. Run the capture recipe, then use the newest scriptlog file.';
}
