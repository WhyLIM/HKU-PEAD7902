# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "pymol-open-source-whl",
# ]
# ///
"""Step 7 - render the two structure views for Figure 6 with PyMOL (headless, OSMesa).

Panel a: the p53 core domain engaged with the DNA duplex (chain B of 1TUP),
grey cartoon, orange DNA, purple Zn2+, red sticks for the two DNA-contact
hotspots (R248, R273).

Panel b: the same domain rotated by about 120 degrees with the DNA hidden; blue
sticks mark the three buried structural hotspots (R175, G245, R282).

Labels are pseudoatoms placed in front of all geometry (along the camera axis)
so that cartoons and DNA cannot hide them, with dashed leader lines back to the
C-alpha atom of each residue.

Outputs (results/pymol/):
  fig6_view_a.png, fig6_view_b.png, fig6_structure.pse
"""
import os

os.environ["PYOPENGL_PLATFORM"] = "osmesa"

import numpy as np
import pymol  # noqa: E402

pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()

from pymol import cmd  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PDB = os.path.join(ROOT, "data", "structure", "1TUP.pdb")
OUT = os.path.join(ROOT, "results", "pymol")
os.makedirs(OUT, exist_ok=True)

cmd.load(PDB, "p53")
if cmd.count_atoms("all") == 0:
    print("ERROR: structure failed to load")
    cmd.quit()
    raise SystemExit(1)

cmd.remove("solvent")
cmd.remove("chain A+C")
cmd.remove("not (chain B or chain E+F or resn ZN)")

cmd.hide("everything")
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("cartoon_transparency", 0)
cmd.set("ray_shadows", 0)
cmd.set("ambient", 0.35)
cmd.set("specular", 0.15)
cmd.set("antialias", 2)
cmd.set("label_size", 13)
cmd.set("label_font_id", 7)
cmd.set("label_color", "black")
cmd.set("label_outline_color", "white")
cmd.set("label_position", (0, 0, 0))
cmd.set("sphere_scale", 0.40)
cmd.set("stick_radius", 0.18)
cmd.set("dash_color", "grey50")
cmd.set("dash_gap", 0.25)
cmd.set("dash_length", 0.18)
cmd.set_color("dna_orange", [0.816, 0.549, 0.376])
cmd.set_color("protein_grey", [0.686, 0.725, 0.765])
cmd.set_color("contact_red", [0.702, 0.149, 0.118])
cmd.set_color("buried_blue", [0.200, 0.361, 0.506])
cmd.set_color("zinc_purple", [0.482, 0.427, 0.553])

cmd.show("cartoon", "chain B")
cmd.color("protein_grey", "chain B")
cmd.show("cartoon", "chain E+F")
cmd.color("dna_orange", "chain E+F")
cmd.show("spheres", "resn ZN")
cmd.color("zinc_purple", "resn ZN")
cmd.show("sticks", "sidechain and chain B and resi 175+245+248+273+282")
cmd.color("contact_red", "chain B and resi 248+273")
cmd.color("buried_blue", "chain B and resi 175+245+282")

coords = []
cmd.iterate_state(1, "chain B and name CA", "coords.append((x,y,z))",
                  space={"coords": coords})
centre = np.mean(np.array(coords), axis=0)


def label_residue(resi, name, right=0.0, up=0.0):
    """Place a label in front of the scene and connect it to the CA atom."""
    view = cmd.get_view()
    right_v, up_v, fwd_v = (np.array(view[0:3]), np.array(view[3:6]),
                            np.array(view[6:9]))
    ca = np.array(cmd.get_atom_coords(f"chain B and resi {resi} and name CA"))
    pos = ca + right_v * right + up_v * up - fwd_v * 8.0
    obj = f"lab_{resi}"
    cmd.pseudoatom(obj, pos=list(pos))
    cmd.label(obj, f'"{name}"')
    cmd.distance(f"lead_{resi}", f"chain B and resi {resi} and name CA", obj)
    cmd.hide("labels", f"lead_{resi}")


cmd.orient("chain B")
cmd.turn("x", 8)
cmd.turn("y", -25)
cmd.zoom("chain B", 4.6)

# panel a: DNA-contact hotspots
label_residue(248, "R248", right=-6, up=-6)
label_residue(273, "R273", right=-6, up=6)
cmd.png(os.path.join(OUT, "fig6_view_a.png"), width=1500, height=1500,
        dpi=300, ray=1)
cmd.save(os.path.join(OUT, "fig6_structure_a.pse"))

# panel b: buried structural hotspots, DNA hidden
cmd.hide("everything", "chain E+F")
cmd.hide("labels")
cmd.delete("lab_*")
cmd.delete("lead_*")
cmd.turn("y", 115)
cmd.turn("x", 6)
label_residue(175, "R175", right=-6, up=-5)
label_residue(245, "G245", right=6, up=-8)
label_residue(282, "R282", right=6, up=7)
cmd.png(os.path.join(OUT, "fig6_view_b.png"), width=1500, height=1500,
        dpi=300, ray=1)
cmd.save(os.path.join(OUT, "fig6_structure_b.pse"))

print("wrote", os.path.join(OUT, "fig6_view_a.png"))
print("wrote", os.path.join(OUT, "fig6_view_b.png"))
cmd.quit()