// History lives in this browser only (IndexedDB).
let opening;

function open() {
  opening ??= new Promise((resolve, reject) => {
    const request = indexedDB.open('foldit-energy', 1);
    request.onupgradeneeded = () => {
      request.result.createObjectStore('steps', { keyPath: 'key' });
      request.result.createObjectStore('deleted'); // so a re-read log does not resurrect them
      request.result.createObjectStore('settings');
    };
    request.onsuccess = () => resolve(request.result);
    request.onerror = () => reject(request.error);
  });
  return opening;
}

async function transaction(names, mode, work) {
  const db = await open();
  return new Promise((resolve, reject) => {
    const tx = db.transaction(names, mode), request = work(tx);
    tx.oncomplete = () => resolve(request?.result);
    tx.onerror = tx.onabort = () => reject(tx.error ?? new Error('Browser storage failed.'));
  });
}

export async function allSteps() {
  const steps = await transaction('steps', 'readonly', tx => tx.objectStore('steps').getAll());
  return steps.sort((a, b) => a.order - b.order);
}

export const getStep = key => transaction('steps', 'readonly', tx => tx.objectStore('steps').get(key));
export const getSetting = key => transaction('settings', 'readonly', tx => tx.objectStore('settings').get(key));
export const setSetting = (key, value) => transaction('settings', 'readwrite', tx => tx.objectStore('settings').put(value, key));

// Adds only steps that are neither present nor deleted; resolves to how many were added.
// With restore (importing a backup), deleted steps come back as well.
export async function addSteps(steps, { restore = false } = {}) {
  let added = 0;
  await transaction(['steps', 'deleted'], 'readwrite', tx => {
    const table = tx.objectStore('steps'), gone = tx.objectStore('deleted');
    for (const step of new Map(steps.map(s => [s.key, s])).values()) {
      if (restore) gone.delete(step.key);
      const deleted = gone.count(step.key);
      deleted.onsuccess = () => {
        if (deleted.result) return;
        const present = table.count(step.key);
        present.onsuccess = () => { if (!present.result) { table.add(step); added++; } };
      };
    }
  });
  return added;
}

// An update must not bring back a step that was deleted while it was being analysed.
export function updateStep(step) {
  return transaction('steps', 'readwrite', tx => {
    const table = tx.objectStore('steps'), present = table.count(step.key);
    present.onsuccess = () => { if (present.result) table.put(step); };
  });
}

export function deleteSteps(steps) {
  return transaction(['steps', 'deleted'], 'readwrite', tx => {
    for (const { key } of steps) {
      tx.objectStore('steps').delete(key);
      tx.objectStore('deleted').put(1, key);
    }
  });
}
