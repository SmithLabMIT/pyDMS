"""
pyDMS.report

Copyright 2025 Brandon C. Tapia

Licensed under the MIT License
"""

import os
import tempfile
import time
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from importlib.resources import files
from PIL import Image
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    Image as PlatypusImage,
    # ListItem,
    # ListFlowable,
    PageBreak,
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.units import inch

from . import visualize as vis


def wait_for_file(path, timeout=20.0):
    """Waits for file to become available when opening to avoid a "permission
        denied" error.

    Args:
        path: The filepath of the file to try to open
        timeout: The number of seconds to continue trying to open before
            failing

    Returns:
        None
    """

    iter1 = False
    start = time.time()
    while True:
        try:
            with open(path, "rb"):
                return
        except PermissionError:
            if iter1 is False:
                print(f"Waiting for {path} to be available...")
                iter1 = True
            if time.time() - start > timeout:
                raise
            time.sleep(0.5)


def get_scaled_image_dimensions(img_path, max_width=7 * inch, max_height=9 * inch, dpi=300):
    with Image.open(img_path) as img:
        width_px, height_px = img.size

        #  1 pt = 1/72 inch
        width_pt = (width_px / dpi) * 72
        height_pt = (height_px / dpi) * 72

        # Determine the uniform scale factor
        scale_w = max_width / width_pt
        scale_h = max_height / height_pt
        scale = min(scale_w, scale_h, 1.0)  # Don't upscale

        desired_width = width_pt * scale
        desired_height = height_pt * scale

    return desired_width, desired_height


def LFER(gas, tmpdir, outliers=False):
    """Plots and saves the LFER linear fits

    report.LFER should be reserved specifically for printing in the report.
    For general plotting of the LFERs, use visualization.LFER

    Args:
        gas: An instance of the Gas class with Gas.LFER populated
        outliers: whether the plot should contain all outlier data

    Returns:
        The name of the file that was saved
    """

    vis.LFER(gas, outliers)
    name = "LFER_outliers.tmp.png" if outliers else "LFER.tmp.png"
    path = os.path.join(tmpdir, name)
    plt.savefig(path, bbox_inches="tight", dpi=300)
    plt.close()
    return path


def LFER_outliers(gas, tmpdir, outliers=True):
    """Plots and saves the LFER linear fits

    report.LFER should be reserved specifically for printing in the report.
    For general plotting of the LFERs, use visualization.LFER

    Args:
        gas: An instance of the Gas class with Gas.LFER populated
        outliers: whether the plot should contain all outlier data

    Returns:
        The name of the file that was saved
    """

    vis.LFER(gas, outliers)
    name = "LFER_outliers.tmp.png" if outliers else "LFER.tmp.png"
    path = os.path.join(tmpdir, name)
    plt.savefig(path, bbox_inches="tight", dpi=300)
    plt.close()
    return path


def histograms(gas, tmpdir):
    """Plots and saves the histograms from the van't Hoff fits

    report.histograms should be reserved specifically for printing in the
        report.
    For general plotting of the histograms, use visualization.histograms

    Args:
        gas: An instance of the Gas class with Gas.vH populated

    Returns:
        The name of the file that was saved
    """

    vis.histograms(gas)
    path = os.path.join(tmpdir, "hist.tmp.png")
    plt.savefig(path, bbox_inches="tight", dpi=300)
    plt.close()
    return path


def isotherms(gas, tmpdir):
    """Plots and saves the sorption isotherms with DMS parameters

    report.isotherms should be reserved specifically for printing in the
        report.
    For general plotting of the isotherms, use visualization.isotherms

    Args:
        gas: An instance of the Gas class with Gas populated

    Returns:
        The name of the file that was saved
    """

    vis.isotherms(gas)
    path = os.path.join(tmpdir, "isotherms.tmp.png")
    plt.savefig(path, bbox_inches="tight", dpi=300)
    plt.close()
    return path


def heat_of_sorption(gas, tmpdir="."):
    """
    *
    """

    vis.heat_of_sorption(gas)
    path = os.path.join(tmpdir, "heat_of_sorption.tmp.png")
    plt.savefig(path, bbox_inches="tight", dpi=300)
    plt.close()
    return path


def footer(canvas, doc):
    """Creates the canvas objects to place in the footer of each page

    Args:
        canvas: For reportlab internal use required for drawing the footer
        doc: For reportlab internal use required for drawing the footer

    Returns:
        None
    """
    width, height = letter

    canvas.setFont("Helvetica", 10)

    page_number = f"Page {doc.page}"
    canvas.drawCentredString(width / 2, 0.5 * inch, page_number)

    current_date = datetime.now().strftime("%Y-%m-%d")
    footer_text = f"Generated on {current_date}"
    canvas.drawCentredString(width / 2, 0.3 * inch, footer_text)

    script_dir = os.path.dirname(os.path.abspath(__file__))

    image_path = files("pyDMS.images").joinpath("lab-logo-transparent-background.png")

    wait_for_file(image_path)

    with Image.open(image_path) as img:
        width_img, height_img = img.size
        aspect_ratio = width_img / height_img

        desired_width = 1.8 * inch
        desired_height = desired_width / aspect_ratio

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            img.save(tmp.name)
            temp_image_path = tmp.name

    canvas.drawImage(
        temp_image_path,
        width - 2.1 * inch,
        0.05 * inch,
        width=desired_width,
        height=desired_height,
        mask="auto",
    )


def generate(gas, report_name):
    """Creates the report PDF

    Args:
        gas: An instance of the gas class with Gas.LFER and Gas.vH populated

    Returns:
        None
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"Generating report {report_name}")

        local_time = time.localtime()
        formatted_datetime = time.strftime("%Y-%m-%d %H:%M:%S %Z", local_time)

        temp = gas.temp

        C_H = gas.CH
        kD = gas.kD
        b = gas.b

        C_H_err = gas.CH_err
        kD_err = gas.kD_err
        b_err = gas.b_err

        LFER_data = gas.LFER
        settings = gas.settings

        # vH_data = gas.vH

        slope_kd, int_kd, slope_b, int_b = LFER_data.fit

        doc = SimpleDocTemplate(report_name, pagesize=letter)
        styles = getSampleStyleSheet()

        centered_style = ParagraphStyle(
            name="CenteredStyle", parent=styles["Normal"], textColor=colors.grey, alignment=1
        )

        elements = []

        elements.append(Paragraph("pyDMS", styles["Title"]))

        elements.append(Paragraph("Reproducible Dual-Mode Sorption Parameters", centered_style))
        elements.append(Spacer(1, 12))

        elements.append(Paragraph(f"Analysis performed at {formatted_datetime}", styles["Normal"]))

        elements.append(Paragraph("DMS Parameters", styles["Heading2"]))

        col_widths = [1.2 * inch, 1.2 * inch, 4 * inch]

        # Define common table style
        table_style = TableStyle(
            [
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
                ("ALIGN", (0, 0), (-1, -1), "LEFT"),
                ("GRID", (0, 0), (-1, -1), 1, colors.black),
                ("BACKGROUND", (0, 1), (-1, -1), colors.white),
            ]
        )

        header_data = [["Temperature", "Parameter", "Result ± Error"]]
        header_table = Table(header_data, colWidths=col_widths, rowHeights=20)
        header_table.setStyle(table_style)
        header_table.hAlign = "LEFT"
        header_table.vAlign = "MIDDLE"
        elements.append(header_table)

        for i, _ in enumerate(kD):

            label = f"{temp[i]} K"

            kD_text = "k<sub>D</sub>"
            C_H_text = "C'<sub>H</sub>"

            kD_paragraph = Paragraph(kD_text, styles["Normal"])
            C_H_paragraph = Paragraph(C_H_text, styles["Normal"])

            data = [
                [label, kD_paragraph, f"{kD[i]} ± {kD_err[i]}"],
                ["", C_H_paragraph, f"{C_H[i]} ± {C_H_err[i]}"],
                ["", "b", f"{b[i]} ± {b_err[i]}"],
            ]

            table = Table(data, colWidths=col_widths, rowHeights=20)

            table.setStyle(table_style)
            table.hAlign = "LEFT"
            header_table.vAlign = "MIDDLE"

            table.setStyle(TableStyle([("SPAN", (0, 0), (0, 2))]))
            table.setStyle(
                TableStyle(
                    [("ALIGN", (0, 0), (0, 2), "LEFT"), ("VALIGN", (0, 0), (0, 2), "MIDDLE")]
                )
            )

            elements.append(table)

        # Sorption Isotherms
        elements.append(PageBreak())

        elements.append(Paragraph("Sorption Isotherms", styles["Heading2"]))
        elements.append(Spacer(1, 12))

        plot_path = isotherms(gas, tmpdir=tmpdir)

        wait_for_file(plot_path)

        desired_width, desired_height = get_scaled_image_dimensions(
            plot_path, max_width=doc.width - 10, max_height=doc.height - 70, dpi=300
        )
        elements.append(PlatypusImage(plot_path, width=desired_width, height=desired_height))

        elements.append(PageBreak())

        elements.append(Paragraph("Linear Free Energy Relationships", styles["Heading2"]))
        elements.append(Spacer(1, 12))

        data = [
            [
                Paragraph("Heat of Henry Sorption (kJ/mol)", styles["Normal"]),
                Paragraph(
                    f"ΔH<sub>D</sub> = {int_kd:.5G} + {slope_kd:.5G} " "ln(k<sub>D,0</sub>)",
                    styles["Normal"],
                ),
            ],
            [
                Paragraph("Heat of Langmuir Sorption (kJ/mol)", styles["Normal"]),
                Paragraph(
                    f"ΔH<sub>b</sub> = {int_b:.5G} + {slope_b:.5G} " "ln(b<sub>0</sub>)",
                    styles["Normal"],
                ),
            ],
        ]

        table = Table(data, colWidths=[3 * inch, 3 * inch])

        table.hAlign = "LEFT"

        table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.white),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
                    # ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                    ("GRID", (0, 0), (-1, -1), 1, colors.black),
                    ("BACKGROUND", (0, 1), (-1, -1), colors.white),
                ]
            )
        )

        elements.append(table)

        elements.append(Spacer(1, 12))

        plot_path = LFER(gas, tmpdir=tmpdir)

        wait_for_file(plot_path)

        with Image.open(plot_path) as img:
            width, height = img.size
            aspect_ratio = width / height

            desired_width = 7 * inch
            desired_height = desired_width / aspect_ratio

        elements.append(PlatypusImage(plot_path, width=desired_width, height=desired_height))

        plot_path = LFER_outliers(gas, tmpdir=tmpdir)

        wait_for_file(plot_path)

        with Image.open(plot_path) as img:
            width, height = img.size
            aspect_ratio = width / height

            desired_width = 7 * inch
            desired_height = desired_width / aspect_ratio

        elements.append(PlatypusImage(plot_path, width=desired_width, height=desired_height))

        elements.append(PageBreak())

        elements.append(Paragraph("van't Hoff Relationships", styles["Heading2"]))
        elements.append(Spacer(1, 12))

        plot_path = histograms(gas, tmpdir=tmpdir)

        wait_for_file(plot_path)

        desired_width, desired_height = get_scaled_image_dimensions(
            plot_path, max_width=doc.width - 10, max_height=doc.height - 70, dpi=300
        )

        elements.append(PlatypusImage(plot_path, width=desired_width, height=desired_height))

        elements.append(PageBreak())

        elements.append(Paragraph("Sorption Energetics", styles["Heading2"]))

        plot_path = heat_of_sorption(gas, tmpdir=tmpdir)

        wait_for_file(plot_path)

        with Image.open(plot_path) as img:
            width, height = img.size
            aspect_ratio = width / height

            desired_width = 7 * inch
            desired_height = desired_width / aspect_ratio

        elements.append(PlatypusImage(plot_path, width=desired_width, height=desired_height))

        data = [
            [
                "",
                Paragraph("Entropic Prefactor (S<sub>i,0</sub>)", styles["Normal"]),
                Paragraph("ΔH<sub>i</sub> (kJ/mol)", styles["Normal"]),
            ],
            [
                "i = Infinite Dilution",
                f"{gas.analysis.deltaH_S_inf[1]:.3e}±{gas.analysis.deltaH_S_inf_err[1]:.3e}",
                f"{gas.analysis.deltaH_S_inf[0]:.6f} ± {gas.analysis.deltaH_S_inf_err[0]:.6f}",
            ],
            [
                "i = Henry",
                f"{gas.analysis.deltaH_D[1]:.3e} ± {gas.analysis.deltaH_D_err[1]:.3e}",
                f"{gas.analysis.deltaH_D[0]:.6f} ± {gas.analysis.deltaH_D_err[0]:.6f}",
            ],
            [
                "i = Langmuir",
                f"{gas.analysis.deltaH_b[1]:.3e} ± {gas.analysis.deltaH_b_err[1]:.3e}",
                f"{gas.analysis.deltaH_b[0]:.6f} ± {gas.analysis.deltaH_b_err[0]:.6f}",
            ],
        ]

        table = Table(data, colWidths=[inch * 1.5, inch * 2.5, inch * 2.5])
        table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.white),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
                    ("ALIGN", (0, 0), (-1, -1), "LEFT"),
                    ("GRID", (0, 0), (-1, -1), 1, colors.black),
                    ("BACKGROUND", (0, 1), (-1, -1), colors.white),
                ]
            )
        )

        elements.append(table)
        elements.append(PageBreak())

        centered_para = ParagraphStyle(name="CenteredCell", parent=styles["Normal"], alignment=1)

        elements.append(Paragraph("Analysis Settings", styles["Heading2"]))

        data = [
            ["Parameter", "Value"],
            [
                Paragraph("ΔH<sub>D</sub> Initial Guess Range", centered_para),
                f"{settings.get('dHD_guess')}",
            ],
            [
                Paragraph("ΔH<sub>b</sub> Initial Guess Range", centered_para),
                f"{settings.get('dHb_guess')}",
            ],
            [
                Paragraph("k<sub>D,0</sub> Initial Guess Range", centered_para),
                f"{settings.get('kD0_guess')}",
            ],
            [
                Paragraph("b<sub>0</sub> Initial Guess Range", centered_para),
                f"{settings.get('b0_guess')}",
            ],
            [
                Paragraph("ΔH<sub>D</sub> Solver Bounds", centered_para),
                f"{settings.get('dHD_bounds')}",
            ],
            [
                Paragraph("ΔH<sub>b</sub> Solver Bounds", centered_para),
                f"{settings.get('dHb_bounds')}",
            ],
            [
                Paragraph("k<sub>D,0</sub> Solver Bounds", centered_para),
                f"{settings.get('kD0_bounds')}",
            ],
            [
                Paragraph("b<sub>0</sub> Solver Bounds", centered_para),
                f"{settings.get('b0_bounds')}",
            ],
        ]

        table = Table(data, colWidths=[inch * 3, inch * 2])
        table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.white),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
                    ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                    ("GRID", (0, 0), (-1, -1), 1, colors.black),
                    ("BACKGROUND", (0, 1), (-1, -1), colors.white),
                ]
            )
        )

        elements.append(table)

        for i, val in enumerate(np.asarray(settings.get("CH_guess"))):

            data = [[Paragraph(f"C'<sub>H</sub> guess #{i+1}", centered_para), val]]

            table = Table(data, colWidths=[inch * 3, inch * 2])
            table.setStyle(
                TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, 0), colors.white),
                        ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
                        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                        ("GRID", (0, 0), (-1, -1), 1, colors.black),
                        ("BACKGROUND", (0, 1), (-1, -1), colors.white),
                    ]
                )
            )

            elements.append(table)

        for i, val in enumerate(np.asarray(settings.get("CH_bounds"))):

            data = [[Paragraph(f"C'<sub>H</sub> bounds #{i+1}", centered_para), val]]

            table = Table(data, colWidths=[inch * 3, inch * 2])
            table.setStyle(
                TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, 0), colors.white),
                        ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
                        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                        ("GRID", (0, 0), (-1, -1), 1, colors.black),
                        ("BACKGROUND", (0, 1), (-1, -1), colors.white),
                    ]
                )
            )

            elements.append(table)

        data = [
            ["Trials", f"{settings.get('trials')}"],
            ["LFER Solver", f"{settings.get('solver_LFER')}"],
            ["maxiter (LFER solver)", f"{settings.get('maxiter_LFER')}"],
            ["van't Hoff Solver", f"{settings.get('solver_vH')}"],
            ["maxiter (van't Hoff solver)", f"{settings.get('maxiter_vH')}"],
            ["ftol (SLSQP)", f"{settings.get('ftol')}"],
            ["xtol (trust-constr)", f"{settings.get('xtol')}"],
            ["gtol (trust-constr)", f"{settings.get('gtol')}"],
        ]

        table = Table(data, colWidths=[inch * 3, inch * 2])
        table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.white),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
                    ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                    ("GRID", (0, 0), (-1, -1), 1, colors.black),
                    ("BACKGROUND", (0, 1), (-1, -1), colors.white),
                ]
            )
        )

        elements.append(table)

        elements.append(PageBreak())

        elements.append(Paragraph("Questions or Comments?", styles["Heading2"]))
        elements.append(
            Paragraph(
                "Please reach out through the pyDMS GitHub page (*INSERT LINK*)", styles["Normal"]
            )
        )

        elements.append(Paragraph("Citation", styles["Heading2"]))
        elements.append(
            Paragraph(
                "If you used pyDMS in any presented or published work please cite:",
                styles["Normal"],
            )
        )
        elements.append(Paragraph("*INSERT PAPER*", styles["Normal"]))

        elements.append(Paragraph("License", styles["Heading2"]))
        elements.append(
            Paragraph(
                "Copyright 2025 Brandon C. Tapia, Jing Ying Yeo, Pablo Dean, Albert X. Wu, Zachary P. Smith."
            )
        )
        elements.append(Spacer(1, 12))
        elements.append(
            Paragraph("pyDMS is covered under the MIT License (https://opensource.org/license/mit)")
        )

        doc.build(elements, onFirstPage=footer, onLaterPages=footer)

    print("pyDMS successful!")
    print("-----------------------END OF PROGRAM------------------------")
