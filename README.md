# Introduction

`pyband` and `pydos` are two python scripts that take in the VASP calculation results (e.g. OUTCAR and PROCAR) and  convert the results to images. It offers a fast and effective way to preview the calcuated results. The image plotting utilizes `matplotlib` package.

# Examples
## `pyband`

When no arguments is given, `pyband` reads in `OUTCAR` (optionally `KPOINTS`) and find the band information within. It then plot the resulting band structure and save it as `band.png`.

```$ pyband```
![band_with_no_args](examples/band1.png | width=100)
 

```
Usage: pyband [options] arg1 arg2

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -f FILENAME, --file=FILENAME
                        location of OUTCAR
  --procar=PROCAR       location of the PROCAR
  -z EFERMI, --zero=EFERMI
                        energy reference of the band plot
  -o BANDIMAGE, --output=BANDIMAGE
                        output image name, "band.png" by default
  -k KPTS, --kpoints=KPTS
                        kpoint path
  -s FIGSIZE, --size=FIGSIZE
                        figure size of the output plot
  -y YLIM               energy range of the band plot
  --lw=LINEWIDTH        linewidth of the band plot
  --dpi=DPI             resolution of the output image
  --occ=OCC             show the occ of selected atoms to each KS orbital
  --occM=OCCMARKER      the marker used in the plot
  --occMs=OCCMARKERSIZE
                        the size of the marker
  --occMc=OCCMARKERCOLOR
                        the color of the marker

Usage: pydos [options] arg1 arg2

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i PROCAR, --input=PROCAR
                        location of the PROCAR
  -p PDOSATOM, --pdos=PDOSATOM
  -l PDOSLABEL, --label=PDOSLABEL
                        label of the pdos
  -z ZERO, --zero=ZERO  energy reference of the band plot
  --sigma=SIGMA         smearing parameter, default 0.02
  -n NEDOS, --nedos=NEDOS
                        number of point in DOS plot
  -o DOSIMAGE, --output=DOSIMAGE
                        output image name, "dos.png" by default
  -s FIGSIZE, --size=FIGSIZE
                        figure size of the output plot
  -x XLIM               x limit of the dos plot
  -y YLIM               energy range of the band plot
  -e EXTRA              extra energy range of the band plot
  --lw=LINEWIDTH        linewidth of the band plot
  --dpi=DPI             resolution of the output image
  --tot                 show total dos
  --notot               not show total dos
  -q, --quiet           not show the resulting image

```
