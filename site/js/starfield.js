// Argus starfield — vanilla canvas. No deps.
//
// Three layers of stars at increasing brightness/size, each with subtle
// twinkle. Respects prefers-reduced-motion.

(() => {
  const canvas = document.getElementById('starfield');
  if (!canvas) return;
  const ctx = canvas.getContext('2d');

  const reduce = window.matchMedia &&
    window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  const DPR = Math.min(window.devicePixelRatio || 1, 2);
  let width = 0, height = 0;

  // star palette — slightly off-white into the cool/warm range
  const palette = [
    'rgba(245, 247, 255, 1)',   // pure starlight
    'rgba(192, 132, 252, 1)',   // purple primary
    'rgba(103, 232, 249, 1)',   // spectral cyan
    'rgba(229, 233, 247, 1)',   // soft white
  ];

  // Three layers — far/mid/near
  const layers = [
    { count: 220, size: [0.4, 1.0], speed: 0.5, twinkle: 0.005 },
    { count: 90,  size: [0.7, 1.6], speed: 1.0, twinkle: 0.012 },
    { count: 35,  size: [1.0, 2.2], speed: 1.6, twinkle: 0.025 },
  ];

  let stars = [];

  function resize() {
    const w = window.innerWidth;
    const h = window.innerHeight;
    width = w; height = h;
    canvas.width  = Math.floor(w * DPR);
    canvas.height = Math.floor(h * DPR);
    canvas.style.width  = w + 'px';
    canvas.style.height = h + 'px';
    ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
    seed();
  }

  function rand(min, max) { return min + Math.random() * (max - min); }

  function seed() {
    stars = [];
    for (const layer of layers) {
      for (let i = 0; i < layer.count; ++i) {
        stars.push({
          x: Math.random() * width,
          y: Math.random() * height,
          r: rand(layer.size[0], layer.size[1]),
          a: rand(0.35, 1.0),
          phase: Math.random() * Math.PI * 2,
          color: palette[Math.floor(Math.random() * palette.length)],
          layer,
        });
      }
    }
  }

  let last = 0;
  function frame(t) {
    const dt = Math.min(40, t - last) || 16;
    last = t;

    ctx.clearRect(0, 0, width, height);

    for (const s of stars) {
      // gentle twinkle via sinusoidal alpha modulation
      s.phase += s.layer.twinkle * dt;
      const tw = 0.55 + 0.45 * Math.sin(s.phase);
      const alpha = s.a * tw;
      ctx.globalAlpha = alpha;
      ctx.fillStyle = s.color;

      // soft glow halo for the brightest layer
      if (s.layer === layers[2]) {
        ctx.beginPath();
        ctx.arc(s.x, s.y, s.r * 3.4, 0, Math.PI * 2);
        const grad = ctx.createRadialGradient(s.x, s.y, 0, s.x, s.y, s.r * 3.4);
        grad.addColorStop(0, s.color.replace('1)', '0.45)'));
        grad.addColorStop(1, s.color.replace('1)', '0)'));
        ctx.fillStyle = grad;
        ctx.fill();
      }

      // crisp star core
      ctx.beginPath();
      ctx.arc(s.x, s.y, s.r, 0, Math.PI * 2);
      ctx.fillStyle = s.color;
      ctx.fill();
    }
    ctx.globalAlpha = 1;

    if (!reduce) requestAnimationFrame(frame);
  }

  // also draw a couple of static "bright" stars with diffraction spikes for
  // that JWST/Hubble look. Kept minimal — too many becomes noisy.
  function drawSpikedStars() {
    const big = [
      { x: width * 0.18, y: height * 0.22, color: 'rgba(192, 132, 252, 1)', r: 1.6 },
      { x: width * 0.82, y: height * 0.16, color: 'rgba(103, 232, 249, 1)', r: 1.4 },
      { x: width * 0.42, y: height * 0.78, color: 'rgba(245, 247, 255, 1)', r: 2.0 },
    ];
    for (const s of big) {
      ctx.save();
      ctx.translate(s.x, s.y);
      ctx.globalAlpha = 0.55;

      // four-point diffraction spike
      const len = 22 + s.r * 10;
      ctx.strokeStyle = s.color;
      ctx.lineWidth = 0.55;
      ctx.beginPath();
      ctx.moveTo(-len, 0); ctx.lineTo(len, 0);
      ctx.moveTo(0, -len); ctx.lineTo(0, len);
      ctx.stroke();

      // bright core
      ctx.globalAlpha = 1;
      ctx.beginPath();
      ctx.arc(0, 0, s.r, 0, Math.PI * 2);
      ctx.fillStyle = s.color;
      ctx.fill();

      ctx.restore();
    }
  }

  window.addEventListener('resize', resize);
  resize();

  if (reduce) {
    // single static frame for reduced motion
    last = performance.now();
    frame(last);
    drawSpikedStars();
  } else {
    requestAnimationFrame((t) => { last = t; frame(t); drawSpikedStars(); });
  }
})();
