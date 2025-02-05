#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create a LaTeX document including SPASSO current run
figures, and some basic information on the Eulerian/Lagrangian diagnostics.

The document includes several sections by default:
    1. Ongoing operations and upcoming stations: to be filled by on-land users 
       after SPASSO figure analysis to suggest up-coming sampling strategy.
    2. Daily Figure analysis: subsection including all the figures produced
       during SPASSO run
Created on Wed Jan  4 09:42:19 2023

@author: lrousselet
"""
import GlobalVars, Library, Functions
import glob, re, os, shutil
from tabulate import tabulate
from pylatex import Document, Section, Subsection, Command, Center, Figure
from pylatex.utils import NoEscape, bold
from pathlib import Path

def fill_document(doc):
    """Add a section, a subsection and some text to the document.
    :param doc: the document
    :type doc: :class:`pylatex.document.Document` instance
    """
    with doc.create(Center()) as centered:
        centered.append('**********************************************')
        centered.append(bold('\nExecutive Summary'))
        centered.append('\n\n Type here your executive summary\n')
        centered.append('**********************************************')
    
    with doc.create(Section('Ongoing operations and upcoming stations')):
        if 'swotCOn' in GlobalVars.Bull:
            nam,lon,lat = GlobalVars.Bull['swotCOn'],GlobalVars.Bull['swotCOlon'],GlobalVars.Bull['swotCOlat']
            crossover={"name":[],"coordinate":[]}
            for i in range(0,len(GlobalVars.Bull['swotCOn'])):
                crossover["name"].append(nam[i])
                crossover["coordinate"].append([lat[i],lon[i]])
            df=Functions.SWOT_passing_time(crossover)
            doc.append('SWOT passing time (UTC) over:\n')
            doc.append(tabulate(df, headers='keys', tablefmt='pipe',showindex=False)+'\n\n')
        doc.append('Type here.')
            
    with doc.create(Section('Daily figures analysis')):
        with doc.create(Subsection('Altimetry, derived currents')):
            allfiles = glob.glob(os.path.join(GlobalVars.Dir['dir_wrk'], '*PHY*.png'))
            files = list(filter(lambda x: not re.search('OW|KE|FTLE|LLADV|ADV|TIMEFROMBATHY', x),allfiles))
            files.sort(key=os.path.getmtime)
            diagfiles = list(set(allfiles)-set(files))
            diagfiles.sort(key=os.path.getmtime)
            if not files:
                doc.append('No figures done.')
            else:
                doc.append('Type here.')
                for f in files:
                    with doc.create(Figure(position='h!')) as fig:
                        fig.add_image(f, width='300px')
                doc.append(Command(command='clearpage'))
        
        with doc.create(Subsection('SST analysis')):
            alfiles = glob.glob(os.path.join(GlobalVars.Dir['dir_wrk'], '*SST*.png'))
            files = list(filter(lambda x: not re.search('ADV', x),alfiles))
            if not files:
                doc.append('No figures done.')
            else:
                doc.append('Type here.')
                for f in files:
                    with doc.create(Figure(position='h!')) as fig:
                        fig.add_image(f, width='300px')
                doc.append(Command(command='clearpage'))
            
        with doc.create(Subsection('Chlorophyll analysis')):
            files = glob.glob(os.path.join(GlobalVars.Dir['dir_wrk'], '*CHL*.png'))
            if not files:
                doc.append('No figures done.')
            else:
                doc.append('Type here.')
                for f in files:
                    with doc.create(Figure(position='h!')) as fig:
                        fig.add_image(f, width='300px')
                doc.append(Command(command='clearpage'))

        with doc.create(Subsection('Eulerian/Lagrangian analysis')):
            if not diagfiles:
                doc.append('No figures done.')
            else:
                txt = "Eulerian diagnostics computed with Copernicus_PHY velocities:\n KE: kinetic energy \n OW: Okubo-Weiss parameter\n\n"
                txt2 = " Lagrangian diagnostics computed by seeding Lagrangian particles every 0.02deg and advected for 30 days backward in time with Copernicus_PHY velocities:\n"\
                    + "FTLE: finite time Lyapunov exponents (convergent fronts detection)\n"\
                        + "LLADV: longitude and latitude advection\n"\
                            +"Retention parameter (based on computing the okubo Weiss parameter along a particle trajectory): Detect trapping structures (colorbar = days water parcels have a positive vorticity)\n"\
                                +"Timefrombathy: Water age since last contact with isobath XXm (precised in figure title)\n\n"\
                                    +"More details available at: https://www.swot-adac.org/resources/swot-adac-products-access/ \n\n"
                doc.append(txt)
                doc.append(txt2)
                for f in diagfiles:
                    with doc.create(Figure(position='h!')) as fig:
                        fig.add_image(f, width='300px')
                doc.append(Command(command='clearpage'))
                        
        with doc.create(Subsection('Wave forecast analysis')):
            files = glob.glob(os.path.join(GlobalVars.Dir['dir_wrk'], '*WAVF*.png'))
            if not files:
                doc.append('No figures done.')
            else:
                doc.append('Type here.')
                for f in files:
                    with doc.create(Figure(position='h!')) as fig:
                        fig.add_image(f, width='300px')
                doc.append(Command(command='clearpage'))
    
    with doc.create(Center()):
        doc.append(bold('Acknowledgments'))
    if GlobalVars.Bull['acknow']:
        with open(os.path.join(GlobalVars.Dir['cruise_path'], GlobalVars.Bull['acknow'])) as f:
            contents = f.read()
        doc.append(contents)
    else:
        doc.append('Type here.')

def find_latex_compiler():
    """Try different LaTeX compilers and return the first available one."""

    compilers = ["pdflatex", "latexmk", "xelatex", "lualatex"]
    for compiler in compilers:
        path = shutil.which(compiler)
        if path:
            print(f"Found LaTeX compiler: {compiler} ({path})")
            return path
    print("No LaTeX compiler found! Please install MiKTeX, TeXLive, or MacTeX.")
    return None

def create(cruise):
    # Document with `\maketitle` command activated
    geometry_options = {"margin": "0.5in"}
    doc = Document(geometry_options=geometry_options)

    doc.preamble.append(Command('title', '['+cruise+']: SPASSO Images Analysis'))
    doc.preamble.append(Command('author', ','.join(GlobalVars.Bull['authors'])))
    doc.preamble.append(Command('date', NoEscape(r'\today')))
    doc.append(NoEscape(r'\maketitle'))
    fill_document(doc)

    try:
        pdf_filename = os.path.join(GlobalVars.Dir['dir_wrk'], f"{cruise}_bulletin{GlobalVars.all_dates['today']}")

        # Automatically detect the best LaTeX compiler
        latex_compiler = find_latex_compiler()
        
        if latex_compiler:
            doc.generate_pdf(pdf_filename, clean_tex=False, compiler=latex_compiler)
            print(f"PDF generated successfully: {pdf_filename}.pdf")
        else:
            print("Could not generate PDF: No LaTeX compiler found.")
    except Exception as e:
        print(f"Error generating PDF: {e}")

    # Screen print and copy
    Library.Done(f"{cruise}_bulletin{GlobalVars.all_dates['today']}.pdf created.")
    Library.Done(f"{cruise}_bulletin{GlobalVars.all_dates['today']}.tex created.")

    try:
        src = Path(pdf_filename + ".pdf")
        dest = Path("../Bulletin") / src.name
        shutil.copy(src, dest)
        print(f"Copied {src} to {dest}")
        Library.Done('Copy in Bulletin/ done.')
    except Exception as e:
        print(f"Error copying file: {e}")
