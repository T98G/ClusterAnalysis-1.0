# ClusterAnalysis Version 1.0

This program performs cluster analysis based on atomic coordinates from a PDB trajectory. It computes cluster sizes, adjacency lists, and neighbor statistics for selected atoms across frames in the trajectory.

---

## 📦 Requirements

- A multi-frame PDB file with `END` statements between frames.
- Fortran compiler (e.g., `gfortran`)
- Coordinates must include `CRYST1` lines to define box dimensions.

---

## 🔧 Compilation

```bash
gfortran -o clusterAnalysis clusterAnalysis.f90
```

---

## ▶️ Usage

```bash
./clusterAnalysis -pdb <pdb_file> -radius <cutoff> -selection <residue_name> [options]
```

### Required Flags:

| Flag         | Description                         |
|--------------|-------------------------------------|
| `-pdb`       | Path to the PDB file                |
| `-radius`    | Cutoff radius for neighbor search   |
| `-selection` | Atom name to select (e.g., `P`, `CA`) |

---

### Optional Flags:

| Flag             | Description                                   | Default                  |
|------------------|-----------------------------------------------|--------------------------|
| `-n`             | Max number of neighbors to consider           | `15`                     |
| `-sizes`         | Output file for main cluster sizes            | `size.xvg`               |
| `-counts`        | Output file for total clusters per frame      | `cluster_counts.xvg`     |
| `-hist`          | Output file for cluster size histogram        | `cluster_hist.xvg`       |
| `-neighbors`     | Output file for adjacency lists               | `neighbors.dat`          |
| `-npd`           | Output file for neighbors probability         | `neighbors_prob.xvg`     |
| `-help` / `-h`   | Show usage information and exit               |                          |

---

## 📤 Output Files

- `neighbors.dat`: Adjacency list per atom per frame
- `cluster_counts.xvg`: Number of clusters in each frame
- `cluster_hist.xvg`: Cluster size distribution (probability)
- `size.xvg`: Main (largest) cluster size per frame
- `neighbors_prob.xvg`: Histogram of number of neighbors

---

## 📈 Visualization

Output files are `.xvg` format, compatible with **XMGrace** for plotting:

```bash
xmgrace cluster_counts.xvg
```

---

## 📌 Notes

- Only atoms matching the given `-selection` will be included in the analysis.
- Box periodicity is accounted for using `CRYST1` dimensions in the PDB file.
- Ensure the PDB file has proper frame separation using `END`.
- Cluster Analysis 2.0 Will be able to read Gromacs Format Trajectories (.xtc format)
---

That's it. Run, analyze, and plot.

Let me know if you want a sample PDB to test.
