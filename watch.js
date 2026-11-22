// Follows the script log. Only Chromium browsers can re-read a file the player picked.
import { getSetting, setSetting } from './store.js';

const INTERVAL = 2000;
export const canFollow = 'showOpenFilePicker' in window;

let handle = null, stamp = '', timer = null;

export function follow(chosen, onText, onState) {
  clearTimeout(timer);
  handle = chosen;
  stamp = '';
  const poll = async () => {
    if (handle !== chosen) return;
    try {
      const file = await chosen.getFile(), now = `${file.lastModified}:${file.size}`;
      if (now !== stamp) {
        await onText(await file.text());
        stamp = now;
      }
      onState(`Following ${chosen.name}`);
    } catch (error) {
      if (error.name === 'NotAllowedError') return onState(`Access to ${chosen.name} expired: press Reconnect`, true);
      onState(`Waiting for ${chosen.name} (${error.message})`); // Foldit may be writing it
    }
    if (handle === chosen) timer = setTimeout(poll, INTERVAL);
  };
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
