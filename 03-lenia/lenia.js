// Lenia engine: continuous cellular automaton on a torus, FFT-based convolution.
(function (root) {
  function FFT(n) {
    const levels = Math.log2(n) | 0;
    if (1 << levels !== n) throw new Error('FFT size must be a power of 2');
    const rev = new Uint32Array(n);
    for (let i = 0; i < n; i++) {
      let r = 0;
      for (let b = 0; b < levels; b++) r = (r << 1) | ((i >> b) & 1);
      rev[i] = r;
    }
    const cos = new Float32Array(n / 2), sin = new Float32Array(n / 2);
    for (let i = 0; i < n / 2; i++) { cos[i] = Math.cos(2 * Math.PI * i / n); sin[i] = Math.sin(2 * Math.PI * i / n); }
    // in-place on re/im with offset & stride
    return function transform(re, im, off, stride, inverse) {
      for (let i = 0; i < n; i++) {
        const j = rev[i];
        if (j > i) {
          const a = off + i * stride, b = off + j * stride;
          let t = re[a]; re[a] = re[b]; re[b] = t;
          t = im[a]; im[a] = im[b]; im[b] = t;
        }
      }
      for (let size = 2; size <= n; size <<= 1) {
        const half = size >> 1, step = n / size;
        for (let i = 0; i < n; i += size) {
          for (let j = i, k = 0; j < i + half; j++, k += step) {
            const wr = cos[k], wi = inverse ? sin[k] : -sin[k];
            const a = off + j * stride, b = off + (j + half) * stride;
            const tr = re[b] * wr - im[b] * wi, ti = re[b] * wi + im[b] * wr;
            re[b] = re[a] - tr; im[b] = im[a] - ti;
            re[a] += tr; im[a] += ti;
          }
        }
      }
      if (inverse) for (let i = 0; i < n; i++) { const a = off + i * stride; re[a] /= n; im[a] /= n; }
    };
  }

  class Lenia {
    constructor(cols, rows, opts = {}) {
      this.cols = cols; this.rows = rows; this.n = cols * rows;
      this.R = opts.R ?? 13;
      this.mu = opts.mu ?? 0.15;
      this.sigma = opts.sigma ?? 0.015;
      this.dt = opts.dt ?? 0.1;
      this.A = new Float32Array(this.n);
      this.U = new Float32Array(this.n); // last neighborhood potential
      this.re = new Float32Array(this.n); this.im = new Float32Array(this.n);
      this.kre = new Float32Array(this.n); this.kim = new Float32Array(this.n);
      this.fftRow = FFT(cols); this.fftCol = FFT(rows);
      this.buildKernel();
    }
    fft2(re, im, inverse) {
      for (let y = 0; y < this.rows; y++) this.fftRow(re, im, y * this.cols, 1, inverse);
      for (let x = 0; x < this.cols; x++) this.fftCol(re, im, x, this.cols, inverse);
    }
    buildKernel() {
      const { cols, rows, R } = this;
      const bell = (x, m, s) => Math.exp(-(((x - m) / s) ** 2) / 2);
      this.kre.fill(0); this.kim.fill(0);
      let sum = 0;
      for (let dy = -R; dy <= R; dy++) for (let dx = -R; dx <= R; dx++) {
        const r = Math.hypot(dx, dy) / R;
        if (r >= 1) continue;
        const v = bell(r, 0.5, 0.15);
        const x = (dx + cols) % cols, y = (dy + rows) % rows;
        this.kre[y * cols + x] += v; sum += v;
      }
      for (let i = 0; i < this.n; i++) this.kre[i] /= sum;
      this.fft2(this.kre, this.kim, false);
    }
    // U = K * A
    potential() {
      const { re, im, kre, kim, n } = this;
      re.set(this.A); im.fill(0);
      this.fft2(re, im, false);
      for (let i = 0; i < n; i++) {
        const ar = re[i], ai = im[i];
        re[i] = ar * kre[i] - ai * kim[i];
        im[i] = ar * kim[i] + ai * kre[i];
      }
      this.fft2(re, im, true);
      this.U.set(re);
    }
    step() {
      this.potential();
      const { A, U, n, mu, sigma, dt } = this;
      const inv2s2 = 1 / (2 * sigma * sigma);
      for (let i = 0; i < n; i++) {
        const d = U[i] - mu;
        const g = 2 * Math.exp(-d * d * inv2s2) - 1;
        let v = A[i] + dt * g;
        A[i] = v < 0 ? 0 : v > 1 ? 1 : v;
      }
    }
    mass() { let s = 0; for (let i = 0; i < this.n; i++) s += this.A[i]; return s; }
    clear() { this.A.fill(0); }
  }

  // smooth value noise in [0,1]
  function valueNoise(cols, rows, scale, rnd = Math.random) {
    const nw = Math.ceil(cols / scale) + 2, nh = Math.ceil(rows / scale) + 2;
    const lat = new Float32Array(nw * nh);
    for (let i = 0; i < lat.length; i++) lat[i] = rnd();
    const out = new Float32Array(cols * rows);
    const sm = t => t * t * (3 - 2 * t);
    for (let y = 0; y < rows; y++) for (let x = 0; x < cols; x++) {
      const fx = x / scale, fy = y / scale;
      const x0 = Math.floor(fx), y0 = Math.floor(fy);
      const tx = sm(fx - x0), ty = sm(fy - y0);
      const a = lat[y0 * nw + x0], b = lat[y0 * nw + x0 + 1];
      const c = lat[(y0 + 1) * nw + x0], d = lat[(y0 + 1) * nw + x0 + 1];
      out[y * cols + x] = (a * (1 - tx) + b * tx) * (1 - ty) + (c * (1 - tx) + d * tx) * ty;
    }
    return out;
  }

  // Seed the field with scattered soft blobs of noise
  function randomSeed(L, density = 0.35) {
    const { cols, rows, R } = L;
    const n1 = valueNoise(cols, rows, R * 1.2), n2 = valueNoise(cols, rows, R * 0.5);
    for (let i = 0; i < L.n; i++) {
      const v = 0.7 * n1[i] + 0.3 * n2[i];
      const t = 1 - density;
      L.A[i] = v > t ? Math.min(1, (v - t) / (1 - t) * 1.6) * (0.6 + 0.4 * Math.random()) : 0;
    }
  }

  const api = { Lenia, valueNoise, randomSeed, FFT };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.LeniaLib = api;
})(typeof self !== 'undefined' ? self : this);
