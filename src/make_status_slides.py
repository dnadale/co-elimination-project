"""Build the lab status slide deck: writeups/status_slides_2026-10-12.pptx"""
from pathlib import Path

from pptx import Presentation
from pptx.util import Inches, Pt

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "results/figures"
prs = Presentation()
prs.slide_width, prs.slide_height = Inches(13.33), Inches(7.5)
blank = prs.slide_layouts[6]


def title(slide, text):
    tb = slide.shapes.add_textbox(Inches(0.5), Inches(0.3), Inches(12.3), Inches(0.9))
    p = tb.text_frame.paragraphs[0]
    p.text = text
    p.font.size, p.font.bold = Pt(30), True


def bullets(slide, items, top=1.4, left=0.7, width=12, size=20):
    tb = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(5.5))
    tf = tb.text_frame
    tf.word_wrap = True
    for i, t in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = "• " + t
        p.font.size = Pt(size)
        p.space_after = Pt(10)


def figure_slide(heading, img, caption, tall=False):
    s = prs.slides.add_slide(blank)
    title(s, heading)
    if tall:
        s.shapes.add_picture(str(img), Inches(0.5), Inches(1.2), height=Inches(5.6))
        tb = s.shapes.add_textbox(Inches(6.5), Inches(2), Inches(6.5), Inches(4))
        tb.text_frame.word_wrap = True
        p = tb.text_frame.paragraphs[0]
        p.text = caption
        p.font.size = Pt(20)
    else:
        s.shapes.add_picture(str(img), Inches(0.5), Inches(1.2), width=Inches(12.3))
        tb = s.shapes.add_textbox(Inches(0.7), Inches(6.5), Inches(12), Inches(0.9))
        tb.text_frame.word_wrap = True
        p = tb.text_frame.paragraphs[0]
        p.text = caption
        p.font.size = Pt(16)


s = prs.slides.add_slide(blank)
title(s, "Phenotype-agnostic gene co-elimination screen")
bullets(s, ["Status report, 2026-10-12", "162 Zoonomia mammals, TOGA gene-loss calls, 19,464 genes",
            "Positive control: left-right asymmetry module (MMP21, DAND5, PKD1L1, CIROP, CIROZ)"], top=2.5, size=24)

s = prs.slides.add_slide(blank)
title(s, "Approach")
bullets(s, ["Genotype-first: find genes with shared loss patterns without knowing the phenotype",
            "Build a gene x species status matrix from TOGA calls (162 species x 19,464 genes)",
            "Score all 25.9 million gene pairs for excess shared loss; keep the top 100",
            "Filter histone and physically linked pairs; test the known LR module as a positive control",
            "Next: phylogenetic correction, artifact checks, reverse functional validation"])

figure_slide("Data quality: matrix and QC", FIG / "1_matrix_qc_overview.png",
             "3,153,168 calls, 0 missing cells, source hashes recorded. Ancestry correction not yet applied.")
figure_slide("Screen output: top 100 candidate pairs", FIG / "2_top100_candidates.png",
             "53 histone-histone, 28 physically linked, 19 unlinked pairs. Unlinked pairs are unvalidated.")
figure_slide("Positive control: the LR module co-lost", FIG / "4_LR_module_fisher.png",
             "All 10 pairs significantly co-lost (BH q < 0.001). P-values ignore shared ancestry (exploratory).")
figure_slide("CIROZ loss is clade-clustered", FIG / "5_CIROZ_losses.png",
             "Lost in all cetaceans and ruminants, 23/24 bats, 17/18 carnivorans; 1/27 primates.")
figure_slide("Module loss on the species tree", ROOT / "results/lr_module/02_tree_with_module_status.png",
             "The five LR genes are lost together in cetaceans and ruminants and intact in most primates and "
             "rodents. Red dot = CIROZ lost. Clade colours assigned by genus.", tall=True)
figure_slide("19 unlinked candidate pairs", FIG / "3_unlinked_pairs.png",
             "Candidates for pathway-coherence follow-up (red = involves CIROZ). Not yet validated.")

s = prs.slides.add_slide(blank)
title(s, "Caveats and next steps")
bullets(s, ["Blind screen alone recovered only CIROZ; DAND5, CIROP, MMP21 are intact in ~70-75% of species "
            "and PKD1L1 has 71 uncertain calls",
            "LR test is a positive control with literature-chosen genes, not a discovery",
            "Check CIROZ/PKD1L1 calls for TOGA over-calling and assembly artifacts",
            "Test whether the 19 unlinked pairs are pathway-coherent",
            "Apply phylogenetic correction, then reverse validation (clinical genetics, model organisms)"],
        size=20)

out = ROOT / "writeups/status_slides_2026-10-12.pptx"
prs.save(out)
print(out)
