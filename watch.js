// Follows the script log. Only Chromium browsers can re-read a file the player picked.
import { getSetting, setSetting } from './store.js';

const INTERVAL = 2000;
export const canFollow = 'showOpenFilePicker' in window;

let handle = null, stamp = '', ticker = null;

// Foldit usually covers the browser, and hidden pages may get one timer tick a minute. Foldit also
// overwrites the log on every run, so a slow poll loses captures. Worker timers keep their pace.
function startTicker(tick) {
  const source = `setInterval(() => postMessage(0), ${INTERVAL});`;
  const worker = new Worker(URL.createObjectURL(new Blob([source], { type: 'text/javascript' })));
  worker.onmessage = tick;
  return worker;
}

export function follow(chosen, onText, onState) {
  ticker?.terminate();
  handle = chosen;
  stamp = '';
  let busy = false;
  const poll = async () => {
    if (handle !== chosen || busy) return;
    busy = true;
    try {
      const file = await chosen.getFile(), now = `${file.lastModified}:${file.size}`;
      if (now !== stamp) {
        await onText(await file.text());
        stamp = now;
      }
      onState(`Connected to ${chosen.name}`);
    } catch (error) {
      if (error.name === 'NotAllowedError') {
        ticker?.terminate();
        return onState(`Access to ${chosen.name} expired: press Reconnect`, true);
      }
      onState(`Waiting for ${chosen.name} (${error.message})`); // Foldit may be writing it
    } finally {
      busy = false;
    }
  };
  ticker = startTicker(poll);
  poll();
}

export async function pickLog() {
  const [chosen] = await window.showOpenFilePicker({
    id: 'foldit-log',
    types: [{ description: 'Foldit script log', accept: { 'text/xml': ['.xml'] } }],
  });
  await setSetting('log', chosen).catch(() => {}); // following works even if it cannot be remembered
  return chosen;
}

export async function rememberedLog() {
  if (!canFollow) return null;
  try {
    const saved = await getSetting('log');
    return saved ? { handle: saved, granted: (await saved.queryPermission({ mode: 'read' })) === 'granted' } : null;
  } catch {
    return null;
  }
}
