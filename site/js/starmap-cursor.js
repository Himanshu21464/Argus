// Argus — site-wide star-map cursor overlay.
//
// Sits over the existing #starfield. The pointer anchors a local
// k-nearest-neighbour graph: stars within reach brighten and connect to
// form transient constellation candidates. Named anchor stars carry
// labels and diffraction spikes. Cursor reticle reports RA / Dec.
//
// Subtle enough to coexist with the static starfield backdrop.
// Honors prefers-reduced-motion.

(() => {
  'use strict';

  if (window.__argusStarmapCursor) return;
  window.__argusStarmapCursor = true;

  const reduce = window.matchMedia &&
    window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  if (reduce) return;

  const canvas = document.createElement('canvas');
  canvas.id = 'starmap-cursor';
  canvas.setAttribute('aria-hidden', 'true');
  canvas.style.cssText = [
    'position:fixed', 'inset:0',
    'width:100%', 'height:100%',
    'pointer-events:none',
    'z-index:2', // above #starfield (z=0/1), below content
    'opacity:0.85',
  ].join(';');
  const mount = () => document.body && document.body.appendChild(canvas);
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', mount, { once: true });
  } else { mount(); }

  const ctx = canvas.getContext('2d', { alpha: true });
  let W = 0, H = 0;
  let stars = [], named = [];

  function buildField() {
    stars = [];
    const targetCount = Math.min(280, Math.floor((W * H) / 5500));
    for (let i = 0; i < targetCount; i++) {
      const mag = Math.pow(Math.random(), 2.8);
      stars.push({
        x: Math.random() * W,
        y: Math.random() * H,
        mag,
        size: 0.3 + mag * 1.8,
        hue: Math.random() < 0.85 ? 200 + Math.random() * 30 : 30 + Math.random() * 30,
        tw: Math.random() * Math.PI * 2,
        twSpeed: 0.015 + Math.random() * 0.035,
      });
    }
    const names = ['POLARIS', 'VEGA', 'ALTAIR', 'DENEB', 'RIGEL', 'BETELGEUSE', 'SIRIUS', 'ARCTURUS', 'SPICA', 'ANTARES'];
    named = names.map((n, i) => {
      const cols = 5, rows = 2;
      const col = i % cols, row = (i / cols) | 0;
      return {
        x: (col + 0.5 + (Math.random() - 0.5) * 0.5) * (W / cols),
        y: (row + 0.5 + (Math.random() - 0.5) * 0.5) * (H / rows),
        name: n,
        size: 2.2 + Math.random() * 1.4,
        hue: 190 + Math.random() * 50,
        tw: Math.random() * Math.PI * 2,
      };
    });
  }

  function resize() {
    W = window.innerWidth; H = window.innerHeight;
    canvas.width = W; canvas.height = H;
    buildField();
  }
  resize();
  window.addEventListener('resize', resize);

  const ptr = { x: -1000, y: -1000 };
  let active = false;
  window.addEventListener('pointermove', (e) => {
    ptr.x = e.clientX; ptr.y = e.clientY; active = true;
  }, { passive: true });

  let raf = 0, paused = false;

  function frame(t) {
    if (paused) return;
    ctx.clearRect(0, 0, W, H);

    const cx = ptr.x, cy = ptr.y;
    const reach = 200, reach2 = reach * reach;

    // Pass 1 — stars + collect nearby
    const nearby = [];
    for (const s of stars) {
      s.tw += s.twSpeed;
      const tw = 0.65 + Math.sin(s.tw) * 0.35;
      let baseA = 0.20 + s.mag * 0.45;
      if (active) {
        const dx = s.x - cx, dy = s.y - cy;
        const d2 = dx * dx + dy * dy;
        if (d2 < reach2) {
          const dist = Math.sqrt(d2);
          const op = 1 - dist / reach;
          baseA += op * 0.4;
          nearby.push({ s, dist, op });
        }
      }
      ctx.fillStyle = `hsla(${s.hue}, 30%, 92%, ${baseA * tw})`;
      ctx.beginPath(); ctx.arc(s.x, s.y, s.size * (0.8 + tw * 0.3), 0, Math.PI * 2); ctx.fill();
    }

    // Pass 2 — proximity edges
    if (nearby.length > 1) {
      const edgeR = 130, edgeR2 = edgeR * edgeR;
      ctx.lineWidth = 0.6;
      for (let i = 0; i < nearby.length; i++) {
        const a = nearby[i];
        for (let j = i + 1; j < nearby.length; j++) {
          const b = nearby[j];
          const dx = b.s.x - a.s.x, dy = b.s.y - a.s.y;
          const d2 = dx * dx + dy * dy;
          if (d2 < edgeR2) {
            const d = Math.sqrt(d2);
            const op = (1 - d / edgeR) * a.op * b.op;
            if (op < 0.06) continue;
            ctx.strokeStyle = `rgba(0, 212, 255, ${op * 0.85})`;
            ctx.beginPath(); ctx.moveTo(a.s.x, a.s.y); ctx.lineTo(b.s.x, b.s.y); ctx.stroke();
          }
        }
      }
    }

    // Pass 3 — crosshair on bright nearby stars + nearest finder
    let nearest = null, nearestD = Infinity;
    for (const n of nearby) {
      if (n.op > 0.5) {
        const r = n.s.size + 3 + n.op * 3;
        ctx.strokeStyle = `rgba(176, 149, 255, ${n.op * 0.45})`;
        ctx.lineWidth = 0.6;
        ctx.beginPath();
        ctx.moveTo(n.s.x - r, n.s.y); ctx.lineTo(n.s.x + r, n.s.y);
        ctx.moveTo(n.s.x, n.s.y - r); ctx.lineTo(n.s.x, n.s.y + r);
        ctx.stroke();
      }
      if (n.dist < nearestD) { nearestD = n.dist; nearest = n.s; }
    }

    // Pass 4 — dashed cursor → nearest star
    if (active && nearest && nearestD < reach) {
      ctx.strokeStyle = 'rgba(0, 212, 255, 0.45)';
      ctx.lineWidth = 0.8;
      ctx.setLineDash([4, 4]);
      ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(nearest.x, nearest.y); ctx.stroke();
      ctx.setLineDash([]);
    }

    // Pass 5 — named anchors (always visible, just dimmer when far)
    ctx.font = '8.5px "JetBrains Mono", ui-monospace, monospace';
    ctx.textBaseline = 'middle';
    for (const n of named) {
      n.tw += 0.025;
      const tw = 0.75 + Math.sin(n.tw) * 0.25;
      let op = 0.5;
      if (active) {
        const dx = n.x - cx, dy = n.y - cy;
        const dist = Math.hypot(dx, dy);
        op = Math.min(1, Math.max(0.4, 1 - dist / 500));
      }
      // diffraction spike
      ctx.strokeStyle = `hsla(${n.hue}, 60%, 90%, ${op * 0.35 * tw})`;
      ctx.lineWidth = 0.5;
      const spike = n.size * (2.4 + tw);
      ctx.beginPath();
      ctx.moveTo(n.x - spike, n.y); ctx.lineTo(n.x + spike, n.y);
      ctx.moveTo(n.x, n.y - spike); ctx.lineTo(n.x, n.y + spike);
      ctx.stroke();
      ctx.fillStyle = `hsla(${n.hue}, 70%, 88%, ${op})`;
      ctx.beginPath(); ctx.arc(n.x, n.y, n.size * (0.9 + tw * 0.2), 0, Math.PI * 2); ctx.fill();
      ctx.fillStyle = `rgba(255, 255, 255, ${op * 0.85})`;
      ctx.beginPath(); ctx.arc(n.x, n.y, n.size * 0.5, 0, Math.PI * 2); ctx.fill();
      ctx.fillStyle = `rgba(240, 245, 255, ${op * 0.45})`;
      ctx.fillText(n.name, n.x + n.size + 7, n.y);
    }

    // Pass 6 — cursor reticle + RA/Dec readout
    if (active) {
      ctx.strokeStyle = 'rgba(0, 212, 255, 0.5)';
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(cx, cy, 18, 0, Math.PI * 2); ctx.stroke();
      ctx.beginPath();
      for (let a = 0; a < 4; a++) {
        const ang = a * Math.PI / 2;
        ctx.moveTo(cx + Math.cos(ang) * 12, cy + Math.sin(ang) * 12);
        ctx.lineTo(cx + Math.cos(ang) * 22, cy + Math.sin(ang) * 22);
      }
      ctx.stroke();
      ctx.fillStyle = 'rgba(0, 212, 255, 0.65)';
      const ra = (cx / W * 24).toFixed(2).padStart(5, '0');
      const dec = ((0.5 - cy / H) * 180).toFixed(2);
      ctx.fillText(`RA  ${ra}h`, cx + 28, cy - 6);
      ctx.fillText(`DEC ${dec}°`, cx + 28, cy + 6);
    }

    raf = requestAnimationFrame(frame);
  }

  document.addEventListener('visibilitychange', () => {
    paused = document.hidden;
    if (!paused) { raf = requestAnimationFrame(frame); }
  });

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', () => { raf = requestAnimationFrame(frame); });
  } else {
    raf = requestAnimationFrame(frame);
  }
})();
