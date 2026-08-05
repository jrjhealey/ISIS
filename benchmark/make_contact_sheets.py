#!/usr/bin/env python3
"""
Tile the four camera aspects of each (molecule, method) render into a single
2x2 contact sheet with a caption, so a whole method can be reviewed at once.

Input:  benchmark/output/renders/{mol}_{method}_{aspect}.png
Output: benchmark/output/sheets/{mol}_{method}.png
"""
import os
from PIL import Image, ImageDraw, ImageFont

RENDERS = 'benchmark/output/renders'
SHEETS = 'benchmark/output/sheets'

MOLECULES = [
    ('spike', 'SARS-CoV-2 spike trimer (6VXX)'),
    ('influenza_np', 'Influenza A nucleoprotein (2IQH)'),
    ('anthrax_lf', 'Anthrax lethal factor (1J7N)'),
]
METHODS = [
    ('bcell_consensus', 'B-cell consensus (linear + conformational vote)'),
    ('bcell_discotope', 'B-cell conformational - DiscoTope'),
    ('tcell_mhc1', 'T-cell MHC-I binding (HLA-A*02:01)'),
    ('innate_glyco', 'N-glycosylation sequon sites'),
    ('structure_protrusion', 'Protrusion index'),
]
ASPECTS = ['front', 'side', 'back', 'top']

CAPTION_H = 84
PAD = 8
THUMB_W = 700


def load_font(size):
    for path in ('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',
                 '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'):
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def main():
    os.makedirs(SHEETS, exist_ok=True)
    made = []

    for mol_tag, mol_name in MOLECULES:
        for meth_tag, meth_name in METHODS:
            paths = [f'{RENDERS}/{mol_tag}_{meth_tag}_{a}.png' for a in ASPECTS]
            if not all(os.path.exists(p) for p in paths):
                print(f'skip {mol_tag}/{meth_tag}: missing renders')
                continue

            thumbs = []
            for p, aspect in zip(paths, ASPECTS):
                im = Image.open(p).convert('RGB')
                scale = THUMB_W / im.width
                im = im.resize((THUMB_W, int(im.height * scale)), Image.LANCZOS)
                d = ImageDraw.Draw(im)
                d.text((14, 10), aspect, fill=(90, 90, 90), font=load_font(26))
                thumbs.append(im)

            tw, th = thumbs[0].size
            sheet = Image.new('RGB',
                              (tw * 2 + PAD * 3, th * 2 + PAD * 3 + CAPTION_H),
                              (255, 255, 255))

            d = ImageDraw.Draw(sheet)
            d.text((PAD + 6, 12), f'{mol_name}', fill=(20, 20, 20), font=load_font(30))
            d.text((PAD + 6, 46), f'{meth_name}', fill=(90, 90, 90), font=load_font(24))

            for i, im in enumerate(thumbs):
                x = PAD + (i % 2) * (tw + PAD)
                y = CAPTION_H + PAD + (i // 2) * (th + PAD)
                sheet.paste(im, (x, y))

            out = f'{SHEETS}/{mol_tag}_{meth_tag}.png'
            sheet.save(out, optimize=True)
            made.append(out)
            print(f'wrote {out}  ({sheet.width}x{sheet.height})')

    print(f'\n{len(made)} contact sheets in {SHEETS}')


if __name__ == '__main__':
    main()
