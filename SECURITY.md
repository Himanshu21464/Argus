# Security policy

Argus is a differentiable astrophysical inference kernel — a C++20
substrate that takes photons (JWST spectra, lensed-image positions,
visibilities) and emits posterior physics. It is a *scientific* kernel,
not a security product, but two of its design commitments — **content-
addressed reproducibility** and **sanitizer-clean memory safety** —
have direct security consequences, and we treat reports against them
seriously.

## Reporting a vulnerability

Preferred channel: **[GitHub Security Advisories](https://github.com/Himanshu21464/Argus/security/advisories/new)**.
A private advisory is the right venue for any finding that affects
deployed runs; it carries built-in coordinated disclosure and credits
the reporter on publication.

Fallback: email **security@argus** (replace with the project security
address once published). PGP-encrypted reports are welcome; the public
key block lives at `/.well-known/pgp-key.asc` on the deployed site.

If you are not sure whether a finding is in scope, send it anyway — we
would rather triage a non-issue than miss one.

## Scope

### In scope

- **Memory safety in the C++ kernel.** Buffer overflows, use-after-
  free, double-free, undefined behaviour reachable from public APIs.
  Tests run under `-fsanitize=address,undefined`; sanitizer reports
  on a public input are treated as vulnerabilities, not bugs.
- **Content-address forgery.** Two distinct retrieval runs that
  produce the same content address; or a tampered input that survives
  the content-address check. The seven retrieval substrate proofs and
  the bit-equal-determinism tests are the load-bearing contract here.
- **Nondeterminism leaks in the kernel.** Wall-clock seeds,
  uninitialized memory, hash randomization, or any other source of
  run-to-run variation on a fixed input + seed. These break the
  reproducibility guarantee and we treat them as security defects.
- **Opacity / lens-model substitution.** A `OpacityKernel` or lens-
  deflection callback that accepts inputs it should reject (e.g. NaN
  cross-sections, negative absorber columns) and corrupts the
  retrieval state silently.
- **Website XSS / CSP bypass.** DOM injection, inline-script smuggling,
  `frame-ancestors` bypass, manifest tampering, or any class of
  finding that would weaken the CSP documented below.
- **Supply chain.** A compromised cdnjs Three.js asset that bypasses
  the SRI pin; a Google Fonts CSS variant that exploits the
  `'unsafe-inline'` style residual.

### Out of scope

- **Numerical accuracy of upstream physics data** (HITRAN, HITEMP,
  ExoMol cross-sections, lens catalogues). These are pulled from
  authoritative archives; their correctness is the upstream's
  responsibility. We *do* defend against malformed files crashing
  the kernel — that's in scope under memory safety.
- **Scientific disagreement.** "Argus retrieved a different abundance
  than paper X" is not a security issue; open a discussion in the
  appropriate channel.
- **Denial of service via legitimate ledger growth, large spectra,
  or expensive HMC chains.** Write a benchmark, not a CVE.
- **Theoretical post-quantum risk.** Argus does not currently sign
  artefacts; once it does (roadmap), this section will narrow.
- **Operational misconfigurations** in someone else's deployment,
  unless the project documentation explicitly recommended the
  misconfiguration.

## Supported versions

| Version | Status              | Window                            |
|---------|---------------------|-----------------------------------|
| 0.7.x   | supported (current) | through v1.0 GA (target Q4 2026)  |
| 0.6.x   | advisories only     | no patches; please upgrade        |
| < 0.6   | end-of-life         | no advisories                     |

## What to expect

| Stage                   | Window                                              |
|-------------------------|-----------------------------------------------------|
| Acknowledgement         | within 48 business hours                            |
| Triage + severity       | within 5 business days                              |
| Fix on `main` + release | critical: 14 days; high: 30; medium: 60             |
| Coordinated disclosure  | 90 days from triage, or 7 days after the patch      |
| Credit                  | with permission, in CHANGELOG.md                    |

We will extend the disclosure window if there is a responsible reason
to (e.g. the fix requires coordinated upstream changes in HITRAN or
Three.js).

## Website security model

The site at `site/index.html` is a single-page brief; the kernel does
not run in the browser. The site is hardened against the usual web
threats (XSS, clickjacking, mixed content, plugin abuse) on the same
defense-in-depth principle as the kernel.

### Layered policy

| Layer                 | Where                | Carries                                                                                  |
|-----------------------|----------------------|------------------------------------------------------------------------------------------|
| HTTP response headers | `site/_headers`      | HSTS, CSP, COOP, CORP, COEP, X-Frame-Options, Permissions-Policy, Accept-CH              |
| HTML `<meta>` tags    | `site/index.html`    | CSP, Referrer-Policy, X-Content-Type-Options, Permissions-Policy (offline / local fallback) |

The two CSP carriers are intentional. HTTP enforces on every response;
`<meta>` is the fallback for local previews and servers that strip
headers. Both are kept in sync.

### Content-Security-Policy (strict)

```
default-src 'self';
script-src 'self' https://cdnjs.cloudflare.com 'inline-speculation-rules';
style-src 'self' 'unsafe-inline' https://fonts.googleapis.com;
font-src 'self' https://fonts.gstatic.com;
img-src 'self' data:;
connect-src 'self';
media-src 'self';
object-src 'none';
manifest-src 'self';
worker-src 'self';
base-uri 'self';
form-action 'self';
frame-ancestors 'none';        ⟵ HTTP only (meta-ignored per CSP spec)
upgrade-insecure-requests
```

What the CSP buys us:

- **No `'unsafe-inline'` on scripts.** The site has zero inline
  executable `<script>` blocks; every JS file (`starfield.js`,
  `transit-3d.js`, `reveal.js`, `starmap-cursor.js`) is served same-
  origin and could be hash-pinned in future via the same mechanism
  Asclepius uses (`asset-manifest.json` + in-browser SHA-256
  verification).
- **`script-src` whitelists `https://cdnjs.cloudflare.com`** *only*
  because Three.js loads from there — and the `<script>` tag pins
  the bytes with `integrity="sha384-..."` (Subresource Integrity).
  A cdnjs compromise that substitutes the file fails SRI and the
  browser refuses to execute it.
- **`object-src 'none'`** eliminates `<embed>` / `<object>` legacy-
  plugin XSS vectors permanently.
- **`base-uri 'self'`** defeats `<base href="https://evil/">`
  injection.
- **`frame-ancestors 'none'`** + `X-Frame-Options: DENY` blocks
  clickjacking.
- **`upgrade-insecure-requests`** auto-upgrades any straggler `http://`
  reference to `https://`.
- **`'unsafe-inline'` on styles remains** — the site has one inline
  `<style>` block (inside `<noscript>`, to keep `reveal`-classed
  sections visible without JS) and one inline `style=` attribute,
  plus Google Fonts CSS varies per User-Agent so SRI is not stable.
  Accepted residual.

### Permissions-Policy

Every capability the brief can't possibly need is denied: camera,
microphone, geolocation, payment, USB, bluetooth, serial,
magnetometer, gyroscope, accelerometer, display-capture, MIDI,
autoplay, picture-in-picture, encrypted-media, idle-detection,
screen-wake-lock, xr-spatial-tracking, interest-cohort (FLoC). Only
`fullscreen` and `web-share` are granted `self`.

### Cross-origin isolation

- `Cross-Origin-Opener-Policy: same-origin` — strict process
  isolation; a tab opened from Argus cannot reach into Argus.
- `Cross-Origin-Resource-Policy: same-origin` — own resources cannot
  be consumed by attacker pages.
- `Cross-Origin-Embedder-Policy: credentialless` — strips cookies
  from cross-origin embeds (no ambient-auth leak); enables high-
  resolution timers cleanly.

### UA-hint suppression

`Accept-CH: ""` — no client hints are sent. A static research brief
does not need to fingerprint the visitor.

### Out of scope (website)

- Network MITM below TLS — out of scope; HSTS preload + COOP / CORP /
  COEP are the maximum we can do from the application.
- A compromised Netlify / Cloudflare egress that swaps both the asset
  and its SRI hash in lock-step — defeats SRI for Three.js and the
  CSP whitelist alike. Signed-manifest deployment (Sigstore-style)
  is on the roadmap.
- Browser zero-days that bypass CSP — out of scope.

## Public-facing security artefacts

- This file (`SECURITY.md`) — the public security policy.
- `site/_headers` — the canonical HTTP security headers, mirrored in
  `site/index.html` meta tags.
- GitHub Security Advisories — the private channel for disclosure.

## Acknowledgments

We list (with permission) every researcher whose disclosure resulted
in a fix in [CHANGELOG.md](CHANGELOG.md). Anonymous reports are
honoured.
