#!/usr/bin/env python3
"""Assemble screenshots into a small looping GIF for the README.

    python3 docs/make_gif.py docs/media/host-scope.gif frame1.jpg frame2.jpg ... [--ms 1800] [--last-ms 3200] [--width 800]

Take the screenshots with docs/demo_blur.js active so vendor names, primer / probe / standard sequences are blurred BEFORE anything is
recorded (see README "화면 미리보기"). Needs Pillow.
"""
import argparse
from pathlib import Path

from PIL import Image


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("frames", nargs="+")
    ap.add_argument("--ms", type=int, default=1800, help="time per frame")
    ap.add_argument("--last-ms", type=int, default=3200, help="time of the last frame before the loop restarts")
    ap.add_argument("--width", type=int, default=800)
    a = ap.parse_args()
    imgs = []
    for f in a.frames:
        im = Image.open(f).convert("RGB")
        if im.width != a.width:
            im = im.resize((a.width, round(im.height * a.width / im.width)), Image.LANCZOS)
        imgs.append(im)
    h = max(im.height for im in imgs)                    # screenshots differ by a few pixels in height: pad them to one size
    padded = []
    for im in imgs:
        if im.height != h:
            canvas = Image.new("RGB", (im.width, h), (255, 255, 255))
            canvas.paste(im, (0, 0))
            im = canvas
        padded.append(im)
    imgs = padded
    pal = [im.quantize(colors=160, method=Image.MEDIANCUT, dither=Image.Dither.NONE) for im in imgs]
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pal[0].save(out, save_all=True, append_images=pal[1:], duration=[a.ms] * (len(pal) - 1) + [a.last_ms], loop=0, optimize=True, disposal=2)
    print(f"{out}: {len(pal)} frames, {out.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()
