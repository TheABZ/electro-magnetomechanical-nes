"""Reconstruct Fig. 8's red polynomial from the user-provided paper PDF.

Optional dependencies: pip install pymupdf pillow
Usage: python scripts/digitize_coupling.py /path/to/paper.pdf
The paper image is not redistributed; only derived numerical points are saved.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path

import fitz
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pdf', type=Path)
    args = parser.parse_args()
    doc = fitz.open(args.pdf)
    page = doc[5]  # Printed page 6, Figure 8.
    images = page.get_images(full=True)
    candidates = [im for im in images if im[2:4] == (1000, 773)]
    if len(candidates) != 1:
        raise ValueError('Expected one 1000x773 Figure 8 image; inspect this PDF manually')
    im = np.asarray(Image.open(io.BytesIO(doc.extract_image(candidates[0][0])['image'])).convert('RGB')).astype(float)
    # Axis intersections inspected in the native embedded image; pixel origin top left.
    left, right, top, bottom = 128., 970., 16., 649.
    points = []
    for x in range(131, 968):
        # Exclude the legend, which lies above the physical curve (theta > 16).
        column = im[180:649, x]
        red = (column[:, 0] > 120) & (column[:, 0] > 1.5*column[:, 1]) & (column[:, 0] > 1.5*column[:, 2])
        rows = np.flatnonzero(red)
        if rows.size:
            row = np.median(rows) + 180
            points.append(((x-left)/(right-left)*0.08-0.04,
                           (bottom-row)/(bottom-top)*20))
    points = np.asarray(points)
    np.savetxt(ROOT/'data/figure8_digitized.csv', points, delimiter=',',
               header='y_m,theta_Vs_per_m', comments='')
    # Normalize displacement before fitting, enforcing the even powers of Eq. 11.
    scale = 0.04
    A = np.column_stack([(points[:, 0]/scale)**i for i in (0, 2, 4, 6, 8)])
    fitted, *_ = np.linalg.lstsq(A, points[:, 1], rcond=None)
    coeffs = np.zeros(9)
    coeffs[::2] = fitted / scale**np.arange(0, 9, 2)
    residual = A@fitted-points[:, 1]
    metadata = dict(source='Cadiou et al. (2022), Figure 8 red approximation; Eq. 11',
                    doi='10.1016/j.jsv.2021.116508',
                    status='reconstructed from raster figure, not author-supplied coefficients',
                    pdf_sha256=hashlib.sha256(args.pdf.read_bytes()).hexdigest(),
                    image_size_px=[1000, 773], axes_px=[left, right, top, bottom],
                    valid_y_m=[-0.04, 0.04], theta_coeffs=coeffs.tolist(),
                    normalized_even_coeffs=fitted.tolist(), scale_m=scale,
                    samples=len(points), rmse_Vs_per_m=float(np.sqrt(np.mean(residual**2))),
                    max_residual_Vs_per_m=float(np.max(np.abs(residual))),
                    uncertainty='Pixel calibration and line thickness are not quantified confidence intervals')
    (ROOT/'data/coupling_fit.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
