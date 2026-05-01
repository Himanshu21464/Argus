// Argus — 3D transit scene.
//
// Three.js WebGL scene: a sun-like star at the centre with a hot
// Jupiter orbiting it. Periodically the planet transits in front of
// the star (from the camera's perspective), at which point we apply a
// brief dimming to the star — the transmission-spectroscopy event the
// kernel exists to invert.

(() => {
  const container = document.getElementById('transit-3d');
  if (!container) return;
  if (typeof window.THREE === 'undefined') return;

  const THREE = window.THREE;
  const reduce = window.matchMedia &&
                 window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  // ─── scene + camera ────────────────────────────────────────────────
  const scene = new THREE.Scene();
  const W0 = container.clientWidth;
  const H0 = container.clientHeight || 480;
  const camera = new THREE.PerspectiveCamera(45, W0 / H0, 0.1, 1000);
  camera.position.set(0, 1.5, 28);
  camera.lookAt(0, 0, 0);

  // ─── renderer ──────────────────────────────────────────────────────
  const renderer = new THREE.WebGLRenderer({
    alpha: true,
    antialias: true,
  });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.setSize(W0, H0);
  container.appendChild(renderer.domElement);
  renderer.domElement.style.display = 'block';

  // ─── star (host) ───────────────────────────────────────────────────
  const starColor   = new THREE.Color(0xffd380);
  const dimmedColor = new THREE.Color(0xffb04e);

  const star = new THREE.Mesh(
    new THREE.SphereGeometry(4.0, 64, 64),
    new THREE.MeshBasicMaterial({ color: starColor.clone() })
  );
  scene.add(star);

  // soft halo (back-side sphere) for that "stellar limb" feel
  const halo = new THREE.Mesh(
    new THREE.SphereGeometry(5.4, 32, 32),
    new THREE.MeshBasicMaterial({
      color: 0xffd380, transparent: true, opacity: 0.12, side: THREE.BackSide
    })
  );
  scene.add(halo);

  // outer corona (very faint, larger)
  const corona = new THREE.Mesh(
    new THREE.SphereGeometry(7.2, 32, 32),
    new THREE.MeshBasicMaterial({
      color: 0xffe080, transparent: true, opacity: 0.04, side: THREE.BackSide
    })
  );
  scene.add(corona);

  // ─── planet (hot Jupiter) ──────────────────────────────────────────
  const planet = new THREE.Mesh(
    new THREE.SphereGeometry(0.95, 64, 64),
    new THREE.MeshStandardMaterial({
      color: 0x1a3a52, roughness: 0.55, metalness: 0.08, emissive: 0x080820,
    })
  );
  scene.add(planet);

  // planet atmosphere — back-side cyan sphere = limb glow
  const atmosphere = new THREE.Mesh(
    new THREE.SphereGeometry(1.15, 64, 64),
    new THREE.MeshBasicMaterial({
      color: 0x00d4ff, transparent: true, opacity: 0.18, side: THREE.BackSide
    })
  );
  scene.add(atmosphere);

  // ─── lights ────────────────────────────────────────────────────────
  const starLight = new THREE.PointLight(0xfff0c0, 2.4, 60, 1.0);
  starLight.position.set(0, 0, 0);
  scene.add(starLight);
  scene.add(new THREE.AmbientLight(0x303048, 0.30));

  // ─── starfield background (pure points) ────────────────────────────
  const bgStarsGeom = new THREE.BufferGeometry();
  const bgStarPositions = [];
  for (let i = 0; i < 800; ++i) {
    const r = 80 + Math.random() * 140;
    const theta = Math.random() * Math.PI * 2;
    const phi = (Math.random() - 0.5) * Math.PI;
    bgStarPositions.push(
      r * Math.cos(theta) * Math.cos(phi),
      r * Math.sin(phi),
      r * Math.sin(theta) * Math.cos(phi)
    );
  }
  bgStarsGeom.setAttribute('position',
    new THREE.Float32BufferAttribute(bgStarPositions, 3));
  const bgStars = new THREE.Points(
    bgStarsGeom,
    new THREE.PointsMaterial({
      color: 0xf0f5ff, size: 0.9, sizeAttenuation: true,
      transparent: true, opacity: 0.7
    })
  );
  scene.add(bgStars);

  // ─── animation ─────────────────────────────────────────────────────
  let t = 0;
  const ORBIT_R = 11;
  const ORBIT_SPEED = 0.0035;

  function frame() {
    t += ORBIT_SPEED;

    // planet on a slightly inclined orbit
    const x = ORBIT_R * Math.cos(t);
    const y = 0.3 * Math.sin(t * 0.7);            // tiny y-wobble
    const z = ORBIT_R * Math.sin(t);
    planet.position.set(x, y, z);
    planet.rotation.y = t * 0.6;
    atmosphere.position.copy(planet.position);

    // detect transit: planet "in front" of star (positive z toward camera)
    // and within the projected stellar disk on the (x, y) plane.
    const projDist = Math.sqrt(x * x + y * y);
    const inTransit = z > 0 && projDist < (4.0 + 0.95);
    const dimAmount = inTransit
      ? Math.max(0, 1.0 - projDist / (4.0 + 0.95))
      : 0;
    star.material.color.lerpColors(starColor, dimmedColor, dimAmount * 0.6);
    halo.material.opacity = 0.12 - 0.04 * dimAmount;

    // gentle starfield rotation
    bgStars.rotation.y = t * 0.05;

    renderer.render(scene, camera);
    if (!reduce) requestAnimationFrame(frame);
  }

  // ─── resize handling ───────────────────────────────────────────────
  function resize() {
    const w = container.clientWidth;
    const h = container.clientHeight || 480;
    renderer.setSize(w, h);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  }
  window.addEventListener('resize', resize);
  resize();

  if (reduce) {
    // single static frame for reduced motion
    frame();
  } else {
    requestAnimationFrame(frame);
  }
})();
