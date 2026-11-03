# Living Letters

A small experiment in three steps, each a single page you can open in a browser: cellular automata that run free for a while, then drift into a word, then let go of it.

**Live demo:** https://tyalor.github.io/living-letters/

| # | Folder | What it is |
|---|--------|------------|
| 01 | [`01-conway/`](01-conway/) | Conway's Game of Life. Canvas renderer, toroidal grid, draw with the mouse, drop in gliders, spaceships, a pulsar, a Gosper gun. |
| 02 | [`02-conway-text/`](02-conway-text/) | The same Life, but after a few seconds of running free the cells are nudged, patch by patch, into a hand-lettered word. It is deliberately imperfect: a tilted marker font, ragged edges, stragglers that never settle, a hold that varies. Type your own word. |
| 03 | [`03-lenia/`](03-lenia/) | [Lenia](https://chakazul.github.io/lenia.html), the continuous generalisation of Life. Real creatures (Orbium, Gyrorbium, Scutium, Paraptera, Helicium, Pyroscutium, Circium) glide around, then the whole field is pulled into the word, holds it while it breathes, and releases it back into the rules. |

No build step, no dependencies. Clone and open any `index.html`.

## Controls

Shared across the pages: **Space** play/pause, **N** single step, **R** random / spawn, **C** clear, **T** replay the text intro. Drag on the canvas to draw.

The Lenia page also has a species picker (which sets the growth parameters μ, σ and the time step), a speed slider and a grid-resolution picker (128, 256 or 512 columns).

## How the Lenia one works

`03-lenia/lenia.js` is a dependency-free engine:

- The field is a `Float32Array` on a torus.
- The neighbourhood potential is a convolution with a ring kernel (radius 13, bell-shaped shell), done with a hand-rolled radix-2 FFT so a 512×256 field steps in roughly 20 ms.
- Growth is the standard Lenia Gaussian, `2·exp(−(u−μ)²/2σ²) − 1`, applied with time step `1/T`.

`03-lenia/creatures.js` holds the creature patterns, decoded from Bert Chan's Lenia catalogue (see credits).

The text intro in both 02 and 03 is the same idea: render the word to an offscreen canvas in a marker font, roughen it with value noise, give every cell a "settle time" from smooth noise plus a mild pull toward the centre, and then, while the automaton keeps running underneath, nudge each cell toward its target once its time has come. Nothing is ever forced all at once.

## Credits

- Conway's Game of Life, John Conway, 1970.
- Lenia and the creature patterns: Bert Wang-Chak Chan, [Chakazul/Lenia](https://github.com/Chakazul/Lenia) (MIT). The patterns in `creatures.js` are decoded from that repository's `animals.json`.

## License

MIT, see [LICENSE](LICENSE).
