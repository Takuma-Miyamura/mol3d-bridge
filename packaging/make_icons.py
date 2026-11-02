"""Draw the Mol3D Bridge icon and write icon.png / icon.ico / icon.icns (run once)."""
import math
from pathlib import Path
from PIL import Image, ImageDraw, ImageFilter

S = 1024
OUT = Path(__file__).resolve().parents[1] / "src" / "mol3d_bridge" / "resources"


def lerp(a, b, t):
    return tuple(int(a[i] + (b[i] - a[i]) * t) for i in range(len(a)))


def sphere(size, color):
    im = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    px = im.load()
    r = size / 2
    for y in range(size):
        for x in range(size):
            dx, dy = (x - r + 0.5) / r, (y - r + 0.5) / r
            d2 = dx * dx + dy * dy
            if d2 > 1:
                continue
            dz = math.sqrt(1 - d2)
            lx, ly, lz = -0.45, -0.55, 0.70
            diff = max(0.0, dx * lx + dy * ly + dz * lz)
            spec = max(0.0, dz * 0.95 + (dx * lx + dy * ly) * 0.3) ** 40
            c = [min(255, int(color[i] * (0.28 + 0.8 * diff) + 255 * 0.55 * spec)) for i in range(3)]
            a = 255 if d2 < 0.97 else int(255 * (1 - d2) / 0.03)
            px[x, y] = (*c, a)
    return im


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    bg = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    grad = Image.new("RGBA", (S, S))
    gd = ImageDraw.Draw(grad)
    for y in range(S):
        gd.line([(0, y), (S, y)], fill=lerp((32, 58, 104), (14, 22, 44), y / S) + (255,))
    mask = Image.new("L", (S, S), 0)
    ImageDraw.Draw(mask).rounded_rectangle([40, 40, S - 40, S - 40], radius=210, fill=255)
    bg.paste(grad, (0, 0), mask)
    d = ImageDraw.Draw(bg)
    # a 2D hexagon sketch (the drawing side) fading behind the 3D model (the 3D side)
    cx, cy, R = 330, 640, 175
    hexa = [(cx + R * math.cos(math.radians(60 * k + 30)), cy + R * math.sin(math.radians(60 * k + 30)))
            for k in range(6)]
    d.line(hexa + [hexa[0]], fill=(150, 190, 255, 150), width=22, joint="curve")
    # bridge arrow
    d.arc([300, 210, 760, 670], 200, 300, fill=(255, 190, 70, 230), width=26)
    d.polygon([(640, 205), (720, 250), (630, 300)], fill=(255, 190, 70, 230))
    # ball-and-stick molecule
    atoms = [((640, 560), 120, (150, 150, 158)), ((470, 790), 95, (235, 40, 40)),
             ((850, 690), 100, (60, 95, 245)), ((700, 330), 70, (235, 235, 235)),
             ((850, 470), 64, (235, 235, 235))]
    bonds = [(0, 1), (0, 2), (0, 3), (2, 4)]
    for a, b in bonds:
        (x1, y1), (x2, y2) = atoms[a][0], atoms[b][0]
        d.line([(x1, y1), (x2, y2)], fill=(90, 95, 110, 255), width=58)
        d.line([(x1 - 6, y1 - 6), (x2 - 6, y2 - 6)], fill=(175, 180, 195, 255), width=24)
    for (x, y), r, col in sorted(atoms, key=lambda t: t[1]):
        sp = sphere(2 * r, col)
        shadow = Image.new("RGBA", sp.size, (0, 0, 0, 0))
        shadow.paste((0, 0, 0, 90), (0, 0), sp.split()[3])
        bg.alpha_composite(shadow.filter(ImageFilter.GaussianBlur(14)), (x - r + 14, y - r + 18))
        bg.alpha_composite(sp, (x - r, y - r))
    bg.save(OUT / "icon.png")
    bg.resize((256, 256), Image.LANCZOS).save(OUT / "icon_256.png")
    bg.save(OUT / "icon.ico", sizes=[(16, 16), (24, 24), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)])
    try:
        bg.save(OUT / "icon.icns")
    except Exception as ex:  # noqa: BLE001
        print("icns not written:", ex)
    print("icons written to", OUT)


if __name__ == "__main__":
    main()
