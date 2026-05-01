// Scroll-triggered reveal animations + animated SVG spectrum drawing.
//
// Adds .in class to .reveal elements when they enter the viewport.
// Plays the spectrum fit-line and data-point fade-in animations once
// the figure is in view.

(() => {
  const reduce = window.matchMedia &&
                 window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  // 1. Scroll reveals on .reveal blocks.
  const targets = document.querySelectorAll('.reveal');
  if ('IntersectionObserver' in window && !reduce) {
    const obs = new IntersectionObserver((entries) => {
      entries.forEach((e) => {
        if (e.isIntersecting) {
          e.target.classList.add('in');
          obs.unobserve(e.target);
        }
      });
    }, { rootMargin: '-80px 0px', threshold: 0.05 });
    targets.forEach((t) => obs.observe(t));
  } else {
    targets.forEach((t) => t.classList.add('in'));
  }

  // 2. Animate the spectrum fit-line: set its stroke-dasharray to its own
  //    length so the CSS animation can draw it from 0 -> full length.
  const fitLine = document.querySelector('.spec-svg .fit-line');
  if (fitLine && fitLine.getTotalLength) {
    const len = fitLine.getTotalLength();
    fitLine.style.setProperty('--fit-len', String(len));
    fitLine.style.strokeDasharray = String(len);
    fitLine.style.strokeDashoffset = String(len);
  }

  // 3. Stagger the spectrum data-point fade-in based on x-position.
  const points = document.querySelectorAll('.spec-svg .data-points > g');
  points.forEach((g, i) => {
    g.style.setProperty('--dp-delay', `${0.6 + i * 0.04}s`);
  });

  // 4. Trigger spectrum animation when the figure is visible.
  const fig = document.querySelector('.hero-spectrum');
  if (fig && 'IntersectionObserver' in window && !reduce) {
    const fobs = new IntersectionObserver((entries) => {
      entries.forEach((e) => {
        if (e.isIntersecting) {
          fig.classList.add('animate');
          fobs.unobserve(fig);
        }
      });
    }, { threshold: 0.15 });
    fobs.observe(fig);
  } else if (fig) {
    fig.classList.add('animate');
  }
})();
