// Runs the physics in a worker so the page stays responsive.
let worker, last = 0;
const waiting = new Map();

function start() {
  // A worker loaded by URL takes its security policy from HTTP headers, which GitHub Pages cannot
  // set. A blob worker inherits the page's policy, so `connect-src 'none'` covers the physics too.
  const source = `import ${JSON.stringify(new URL('./worker.js', import.meta.url).href)};`;
  worker = new Worker(URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), { type: 'module' });
  worker.onmessage = ({ data: { id, result, error, progress } }) => {
    const job = waiting.get(id);
    if (!job) return;
    if (progress !== undefined) return job.onProgress?.(progress);
    waiting.delete(id);
    if (error !== undefined) job.reject(new Error(error));
    else job.resolve(result);
  };
  worker.onerror = event => {
    for (const job of waiting.values()) job.reject(new Error(event.message || 'The analysis engine failed to load or ran out of memory.'));
    waiting.clear();
    worker.terminate();
    worker = null;
  };
}

export function run(job, data, onProgress) {
  if (!worker) start();
  return new Promise((resolve, reject) => {
    const id = ++last;
    waiting.set(id, { resolve, reject, onProgress });
    worker.postMessage({ id, job, data });
  });
}
