// Argus starfield — dim background field of stars.
//
// Intentionally restrained: most of the page's astronomy comes from
// the inline SVG (transmission spectrum, atmospheric cross-section,
// constellation chart). The starfield exists only to give the deep-void
// background a sense of scale and depth without competing with the data
// graphics in front.
//
// Honors prefers-reduced-motion.

(() => {
  // create canvas dynamically so the markup stays clean
  const canvas = document.createElement('canvas');
  canvas.id = 'starfield';
  canvas.setAttribute('aria-hidden', 'true');
  document.body.appendChild(canvas);
  const ctx = canvas.getContext('2d');

  const reduce = window.matchMedia &&
    window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  const DPR = Math.min(window.devicePixelRatio || 1, 2);
  let width = 0, height = 0;

  // dim, cool palette — nothing flashy
  const palette = [
    'rgba(240, 245, 255, 1)',  // starlight
    'rgba(184, 194, 212, 1)',  // soft white
    'rgba(103, 232, 249, 1)',  // h2o cyan, rare
    'rgba(176, 149, 255, 1)',  // continuum violet, rare
  ];

  // sparser than typical to keep the data graphics in front readable
  const layers = [
    { count: 140, size: [0.25, 0.85], twinkle: 0.004 },
    { count: 60,  size: [0.55, 1.25], twinkle: 0.010 },
  ];

  let stars = [];

  function resize() {
    const w = window.innerWidth;
    // height = full document height so the field doesn't end at fold
    const h = Math.max(window.innerHeight, document.documentElement.scrollHeight);
    width = w; height = h;
    canvas.width = Math.floor(w * DPR);
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
      const dens = (width * height) / (1440 * 900);
      const n = Math.max(40, Math.floor(layer.count * dens));
      for (let i = 0; i < n; ++i) {
        const colorIdx = Math.random() < 0.85 ? 0 : (1 + Math.floor(Math.random() * 3));
        stars.push({
          x: Math.random() * width,
          y: Math.random() * height,
          r: rand(layer.size[0], layer.size[1]),
          baseAlpha: rand(0.22, 0.65),
          phase: Math.random() * Math.PI * 2,
          color: palette[Math.min(colorIdx, palette.length - 1)],
          twinkle: layer.twinkle,
        });
      }
    }
  }

  function drawSpike(x, y, color, len) {
    ctx.save();
    ctx.translate(x, y);
    ctx.globalAlpha = 0.4;
    ctx.strokeStyle = color;
    ctx.lineWidth = 0.5;
    ctx.beginPath();
    ctx.moveTo(-len, 0); ctx.lineTo(len, 0);
    ctx.moveTo(0, -len); ctx.lineTo(0, len);
    ctx.stroke();
    ctx.restore();
  }

  let last = 0;
  function frame(t) {
    const dt = Math.min(40, t - last) || 16;
    last = t;
    ctx.clearRect(0, 0, width, height);

    for (const s of stars) {
      s.phase += s.twinkle * dt;
      const tw = 0.6 + 0.4 * Math.sin(s.phase);
      ctx.globalAlpha = s.baseAlpha * tw;
      ctx.fillStyle = s.color;
      ctx.beginPath();
      ctx.arc(s.x, s.y, s.r, 0, Math.PI * 2);
      ctx.fill();
    }
    ctx.globalAlpha = 1;

    if (!reduce) requestAnimationFrame(frame);
  }

  function drawSpikes() {
    // a few brighter "JWST-style" stars with diffraction spikes,
    // placed deterministically so they don't shift on resize
    const big = [
      { x: width * 0.08,  y: height * 0.08, color: 'rgba(176, 149, 255, 1)', len: 22 },
      { x: width * 0.93,  y: height * 0.04, color: 'rgba(103, 232, 249, 1)', len: 18 },
      { x: width * 0.36,  y: height * 0.30, color: 'rgba(240, 245, 255, 1)', len: 16 },
      { x: width * 0.72,  y: height * 0.55, color: 'rgba(176, 149, 255, 1)', len: 20 },
    ];
    for (const s of big) {
      drawSpike(s.x, s.y, s.color, s.len);
      ctx.globalAlpha = 1;
      ctx.beginPath();
      ctx.arc(s.x, s.y, 1.4, 0, Math.PI * 2);
      ctx.fillStyle = s.color;
      ctx.fill();
    }
  }

  window.addEventListener('resize', resize);

  // re-seed when the layout settles (images / fonts loaded change page height)
  let scheduled = false;
  function scheduleResize() {
    if (scheduled) return;
    scheduled = true;
    requestAnimationFrame(() => { scheduled = false; resize(); drawSpikes(); });
  }
  window.addEventListener('load', scheduleResize);
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(scheduleResize);

  resize();

  if (reduce) {
    last = performance.now();
    frame(last);
    drawSpikes();
  } else {
    requestAnimationFrame((t) => { last = t; frame(t); drawSpikes(); });
  }
})();
