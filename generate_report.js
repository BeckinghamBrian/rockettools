/**
 * generate_report.js
 * ==================
 * Builds an Engine Analysis Report (.docx) from the JSON data file written by
 * Engine_Class.generate_report().
 *
 * Usage (called automatically by Python):
 *   node generate_report.js <path/to/report_data.json>
 *
 * Output:
 *   <out_dir>/<report_name>.docx
 *
 * Requires: the 'docx' npm package (globally installed or in node_modules).
 *
 * DOCUMENT STRUCTURE
 * ------------------
 *  Cover line (engine name + date)
 *  ─────────────────────────────────
 *  SECTION 1 · INPUTS
 *    Two-column table: Parameter | Value
 *    Followed by a TOML block — copy-paste directly into a new input file.
 *  ─────────────────────────────────
 *  SECTION 2 · CHAMBER INFO & DIMENSIONS
 *    Two-column table: Parameter | Value
 *  ─────────────────────────────────
 *  SECTION 3 · PLOTS
 *    Two plots per page, each 9 × 4 in (fits US Letter with 1-in margins).
 */

'use strict';

const fs   = require('fs');
const path = require('path');

const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  ImageRun, Header, Footer, PageNumber, PageBreak,
  AlignmentType, BorderStyle, WidthType, ShadingType, VerticalAlign,
  PageOrientation,
} = require('docx');

// ── Design constants ──────────────────────────────────────────────────────────
const BLUE        = '2E5FA3';
const LIGHT_BLUE  = 'D6E4F7';
const GREY_ROW    = 'F5F5F5';

// US Letter, 1-inch margins, all values in DXA (1440 DXA = 1 inch)
const PAGE_W    = 12240;
const PAGE_H    = 15840;
const MARGIN    = 1440;
const CONTENT_W = PAGE_W - 2 * MARGIN;   // 9360 DXA  = 6.5 in

// Landscape plot pages: docx-js convention — pass portrait dimensions and set
// orientation: LANDSCAPE; docx-js swaps width/height in the XML internally.
// Content width on landscape page = long edge − 2 margins = 15840 − 2880 = 12960 DXA
const LAND_CONTENT_W = PAGE_H - 2 * MARGIN;   // 12960 DXA  = 9 in
// Content height on landscape page = short edge − 2 margins = 12240 − 2880 = 9360 DXA
const LAND_CONTENT_H = PAGE_W - 2 * MARGIN;   // 9360 DXA  = 6.5 in

// Image dimensions for landscape pages: fill the full landscape content width.
// Each plot PNG was saved at 9 × 4 in by Python.  On landscape we have 9 in
// of width, so we keep width the same and the aspect ratio handles height.
// At 96 dpi: 9 in = 864 px wide, 4 in = 384 px tall — fills half the page.
const IMG_W_PX = 864;   // plot image width in docx pixel units (9 in equivalent)
const IMG_H_PX = 384;   // plot image height in docx pixel units (4 in equivalent)

// ── Utility builders ──────────────────────────────────────────────────────────

const cellBorder = { style: BorderStyle.SINGLE, size: 1, color: 'CCCCCC' };
const allBorders = {
  top: cellBorder, bottom: cellBorder,
  left: cellBorder, right: cellBorder,
};

function heading(text, size = 32, spacing = 360) {
  return new Paragraph({
    spacing: { before: spacing, after: 120 },
    border:  size >= 32
      ? { bottom: { style: BorderStyle.SINGLE, size: 6, color: BLUE, space: 1 } }
      : {},
    children: [new TextRun({ text, bold: true, size, color: BLUE, font: 'Arial' })],
  });
}

function para(text, opts = {}) {
  return new Paragraph({
    spacing: { before: 60, after: 60 },
    children: [new TextRun({ text, font: 'Arial', size: 20, ...opts })],
  });
}

function mono(text) {
  // Monospace paragraph used for the TOML block
  return new Paragraph({
    spacing: { before: 0, after: 0 },
    children: [new TextRun({ text, font: 'Courier New', size: 18 })],
  });
}

function spacer(before = 120) {
  return new Paragraph({ spacing: { before, after: 0 }, children: [] });
}

function pageBreak() {
  return new Paragraph({ children: [new PageBreak()] });
}

function makeCell(text, width, { fill = 'FFFFFF', bold = false } = {}) {
  return new TableCell({
    borders:       allBorders,
    width:         { size: width, type: WidthType.DXA },
    shading:       { fill, type: ShadingType.CLEAR },
    margins:       { top: 80, bottom: 80, left: 120, right: 120 },
    verticalAlign: VerticalAlign.CENTER,
    children: [new Paragraph({
      children: [new TextRun({ text: String(text), font: 'Arial', size: 18, bold })],
    })],
  });
}

/** Two-column table from [[label, value], ...] */
function twoColTable(rows) {
  const colW1 = Math.round(CONTENT_W * 0.58);
  const colW2 = CONTENT_W - colW1;

  const headerRow = new TableRow({
    tableHeader: true,
    children: [
      makeCell('Parameter', colW1, { fill: LIGHT_BLUE, bold: true }),
      makeCell('Value',     colW2, { fill: LIGHT_BLUE, bold: true }),
    ],
  });

  const dataRows = rows.map(([label, value], idx) =>
    new TableRow({
      children: [
        makeCell(label, colW1, { fill: idx % 2 === 0 ? 'FFFFFF' : GREY_ROW }),
        makeCell(value, colW2, { fill: idx % 2 === 0 ? 'FFFFFF' : GREY_ROW }),
      ],
    })
  );

  return new Table({
    width:        { size: CONTENT_W, type: WidthType.DXA },
    columnWidths: [colW1, colW2],
    rows:         [headerRow, ...dataRows],
  });
}

// ── TOML serialiser ───────────────────────────────────────────────────────────
// Produces a syntactically valid TOML file from the inputs object so the reader
// can copy the entire block and use it as a new input file without any edits.

function toToml(inputs) {
  const lines = [
    '# Engine input file',
    '# Generated by Engine_Class — copy this block into a new .toml file',
    '',
    '# ── Required ──────────────────────────────────────────────────────────',
  ];

  const REQUIRED  = ['OF', 'Pc_psi', 'thrust_lbf', 'height_of_optimization', 'burnTime'];
  const LISTS     = ['newMono', 'newFuel', 'newOx'];
  const SECTIONS  = {
    '# ── Design ────────────────────────────────────────────────────────────':
      ['engine_name', 'num_characteristics', 'Lstar', 'alpha', 'beta', 'OF_low', 'OF_high'],
    '# ── Propellant ─────────────────────────────────────────────────────────':
      ['currFuel', 'currOx', 'oxCooled', 'monoMode',
       'newMonoBool', 'newMonoName', 'newMono',
       'newFuelBool', 'newFuelName', 'newFuel',
       'newOxBool',  'newOxName',  'newOx'],
    '# ── Cooling ────────────────────────────────────────────────────────────':
      ['coolantTempStart', 'Kwall', 'Kc', 'chamberWallThickness',
       'channelHeight', 'channelWallThickness', 'numChannels',
       'temp_step', 'coefThermEx', 'youngMod'],
    '# ── Plots ──────────────────────────────────────────────────────────────':
      ['plot_engine_contour', 'plot_nozzle_contour', 'plot_isp_vs_of',
       'plot_mach', 'plot_heat_flux', 'plot_wall_temp', 'plot_thermal_stress',
       'plot_coolant_velocity', 'plot_hx_coefficient', 'plot_coolant_temp',
       'plot_channel_width'],
    '# ── Report ─────────────────────────────────────────────────────────────':
      ['report_output_dir', 'report_dpi'],
  };

  function serVal(val) {
    if (typeof val === 'boolean') return String(val);
    if (typeof val === 'number')  return String(val);
    if (typeof val === 'string')  return `"${val}"`;
    if (Array.isArray(val)) {
      // list-of-lists
      const inner = val.map(row =>
        '    [' + row.map(v => `"${v}"`).join(', ') + ']'
      ).join(',\n');
      return `[\n${inner},\n]`;
    }
    return JSON.stringify(val);
  }

  // Required block
  for (const k of REQUIRED) {
    if (k in inputs) lines.push(`${k} = ${serVal(inputs[k])}`);
  }

  // Sectioned blocks
  for (const [sectionComment, keys] of Object.entries(SECTIONS)) {
    lines.push('');
    lines.push(sectionComment);
    for (const k of keys) {
      if (k in inputs) lines.push(`${k} = ${serVal(inputs[k])}`);
    }
  }

  return lines.join('\n');
}

// ── Label maps ────────────────────────────────────────────────────────────────

const INPUT_LABELS = {
  OF:                    'Mixture Ratio (O/F)',
  Pc_psi:                'Chamber Pressure (psi)',
  thrust_lbf:            'Thrust (lbf)',
  height_of_optimization:'Optimisation Altitude (m)',
  burnTime:              'Burn Time (s)',
  engine_name:           'Engine Name',
  num_characteristics:   'MOC Characteristics',
  Lstar:                 'L* Characteristic Length (m)',
  alpha:                 'Converging Half-angle α (°)',
  beta:                  'Diverging Initial Half-angle β (°)',
  OF_low:                'OF Sweep — Low',
  OF_high:               'OF Sweep — High',
  currFuel:              'Fuel (CEA name)',
  currOx:                'Oxidiser (CEA name)',
  oxCooled:              'Oxidiser-Cooled?',
  coolantTempStart:      'Coolant Inlet Temp (K)',
  Kwall:                 'Wall Thermal Conductivity (W/m·K)',
  Kc:                    'Coolant Thermal Conductivity (W/m·K)',
  chamberWallThickness:  'Wall Thickness (mm)',
  channelHeight:         'Channel Depth (mm)',
  channelWallThickness:  'Channel Wall Thickness (mm)',
  numChannels:           'Number of Cooling Channels',
  temp_step:             'Wall Temp Convergence Step (K)',
  coefThermEx:           'CTE (1/K)',
  youngMod:              'Young\'s Modulus (GPa)',
  monoMode:              'Monoprop Mode',
  newMonoBool:           'Custom Monoprop Enabled',
  newMonoName:           'Custom Monoprop Name',
  newFuelBool:           'Custom Fuel Enabled',
  newFuelName:           'Custom Fuel Name',
  newOxBool:             'Custom Oxidiser Enabled',
  newOxName:             'Custom Oxidiser Name',
  report_output_dir:     'Report Output Directory',
  report_dpi:            'Plot DPI',
  plot_engine_contour:   'Plot — Engine Contour',
  plot_nozzle_contour:   'Plot — Nozzle Contour',
  plot_isp_vs_of:        'Plot — Isp vs OF',
  plot_mach:             'Plot — Mach Number',
  plot_heat_flux:        'Plot — Heat Flux',
  plot_wall_temp:        'Plot — Wall Temperature',
  plot_thermal_stress:   'Plot — Thermal Stress',
  plot_coolant_velocity: 'Plot — Coolant Velocity',
  plot_hx_coefficient:   'Plot — HX Coefficient',
  plot_coolant_temp:     'Plot — Coolant Temperature',
  plot_channel_width:    'Plot — Channel Width',
};

// List-valued keys shown only in TOML, not in the table
const SKIP_IN_TABLE = new Set(['newMono', 'newFuel', 'newOx']);

const CHAMBER_LABELS = {
  Ae_over_At:    'Ae / At  (area ratio)',
  Pc_kPa:        'Chamber Pressure, Pc (kPa)',
  Tc_K:          'Chamber Temperature, Tc (K)',
  thrust_N:      'Thrust (N)',
  Pt_kPa:        'Throat Pressure, Pt (kPa)',
  Tt_K:          'Throat Temperature, Tt (K)',
  Me:            'Exit Mach Number',
  cea_Me:        'CEA Exit Mach Number',
  Po_kPa:        'Stagnation Pressure, Po (kPa)',
  To_K:          'Stagnation Temperature, To (K)',
  Ueq_m_s:       'Exit Velocity, Ueq (m/s)',
  mdot_kg_s:     'Total Mass Flow Rate (kg/s)',
  mdotFuel_kg_s: 'Fuel Mass Flow Rate (kg/s)',
  mdotOx_kg_s:   'Oxidiser Mass Flow Rate (kg/s)',
  Isp_s:         'Specific Impulse, Isp (s)',
  cea_Isp_s:     'CEA Specific Impulse (s)',
  Cstar_m_s:     'Characteristic Velocity, C* (m/s)',
  cp_J_kgK:      'Specific Heat, Cp (J/kg·K)',
  VdotFuel_L_s:  'Fuel Volume Flow Rate (L/s)',
  VdotOx_L_s:    'Oxidiser Volume Flow Rate (L/s)',
  volFuel_L:     'Fuel Tank Volume Required (L)',
  volOx_L:       'Oxidiser Tank Volume Required (L)',
  massFuel_kg:   'Fuel Mass (kg)',
  massOx_kg:     'Oxidiser Mass (kg)',
  massProp_kg:   'Total Propellant Mass (kg)',
  Vc_m3:         'Chamber Volume (m³)',
  Lc_m:          'Chamber Length (m)',
  Dc_m:          'Chamber Diameter (m)',
  At_m2:         'Throat Area (m²)',
  Dt_m:          'Throat Diameter (m)',
  Ae_m2:         'Exit Area (m²)',
  De_m:          'Exit Diameter (m)',
  Rc_mm:         'Chamber Radius, Rc (mm)',
  Rc1_mm:        'Converging Curve Radius 1, Rc1 (mm)',
  Rc2_mm:        'Converging Curve Radius 2, Rc2 (mm)',
  Rt_mm:         'Throat Radius, Rt (mm)',
  Rc3_mm:        'Diverging Curve Radius, Rc3 (mm)',
  Re_mm:         'Exit Radius, Re (mm)',
};

// ── Main ──────────────────────────────────────────────────────────────────────

async function main() {
  const jsonPath = process.argv[2];
  if (!jsonPath) {
    console.error('Usage: node generate_report.js <report_data.json>');
    process.exit(1);
  }

  const data = JSON.parse(fs.readFileSync(jsonPath, 'utf8'));
  const { engine_name, report_name, out_dir, inputs, chamber, plots } = data;

  const children = [];

  // ── Cover ─────────────────────────────────────────────────────────────────
  children.push(new Paragraph({
    spacing: { before: 0, after: 120 },
    children: [new TextRun({
      text: `Engine Analysis Report`,
      bold: true, size: 48, color: BLUE, font: 'Arial',
    })],
  }));
  children.push(new Paragraph({
    spacing: { before: 0, after: 60 },
    children: [new TextRun({
      text: engine_name || 'Unnamed Engine',
      bold: true, size: 32, color: '444444', font: 'Arial',
    })],
  }));
  children.push(para(`Generated: ${new Date().toLocaleString()}`, { color: '888888', size: 18 }));
  children.push(spacer(360));

  // ══════════════════════════════════════════════════════════════════════════
  // SECTION 1 · INPUTS
  // ══════════════════════════════════════════════════════════════════════════
  children.push(heading('1  ·  Inputs', 32, 120));
  children.push(spacer(60));

  // 1a — scalar inputs table
  const inputRows = Object.entries(inputs)
    .filter(([k]) => !SKIP_IN_TABLE.has(k))
    .map(([k, v]) => [INPUT_LABELS[k] || k, String(v)]);
  children.push(twoColTable(inputRows));
  children.push(spacer(300));

  // 1b — TOML block (copy-pasteable new input file)
  children.push(heading('TOML Input Block  (copy → new .toml file)', 24, 120));
  children.push(para(
    'The block below is a syntactically valid TOML file. Copy it into a new ' +
    '.toml file and pass it to Engine.from_file() to reproduce or modify this run.',
    { color: '555555', size: 18 }
  ));
  children.push(spacer(60));
  const tomlLines = toToml(inputs).split('\n');
  tomlLines.forEach(line => children.push(mono(line)));
  children.push(spacer(360));

  // ══════════════════════════════════════════════════════════════════════════
  // SECTION 2 · CHAMBER INFO & DIMENSIONS
  // ══════════════════════════════════════════════════════════════════════════
  children.push(pageBreak());
  children.push(heading('2  ·  Chamber Info & Dimensions', 32, 60));
  children.push(spacer(60));

  const chamberRows = Object.entries(chamber)
    .map(([k, v]) => [CHAMBER_LABELS[k] || k, String(v)]);
  children.push(twoColTable(chamberRows));
  children.push(spacer(360));

  // ══════════════════════════════════════════════════════════════════════════
  // SECTION 3 · PLOTS  (landscape pages, two plots per page)
  //
  // docx supports multiple sections with different page orientations in a
  // single document.  Portrait sections (cover, inputs, chamber) stay as-is.
  // Each PAIR of plots gets its own landscape section.  This means:
  //   - Plots 0+1  → landscape section 1
  //   - Plots 2+3  → landscape section 2  … etc.
  //
  // docx-js landscape convention: pass portrait dimensions (PAGE_W, PAGE_H)
  // with orientation: PageOrientation.LANDSCAPE — docx-js swaps the axes in
  // the XML so the long edge becomes the width on screen and in print.
  //
  // The "no plots" case adds a single note paragraph to the portrait section
  // and does not create any landscape sections.
  // ══════════════════════════════════════════════════════════════════════════

  // ── Shared header/footer builders (reused across all sections) ────────────
  function makeHeader() {
    return new Header({
      children: [new Paragraph({
        border: { bottom: { style: BorderStyle.SINGLE, size: 4, color: 'CCCCCC', space: 1 } },
        children: [new TextRun({
          text: `${engine_name || 'Engine'}  —  Analysis Report`,
          font: 'Arial', size: 18, color: '888888',
        })],
      })],
    });
  }

  function makeFooter() {
    return new Footer({
      children: [new Paragraph({
        alignment: AlignmentType.RIGHT,
        children: [
          new TextRun({ text: 'Page ', font: 'Arial', size: 18, color: '888888' }),
          new TextRun({ children: [PageNumber.CURRENT], font: 'Arial', size: 18, color: '888888' }),
          new TextRun({ text: ' of ', font: 'Arial', size: 18, color: '888888' }),
          new TextRun({ children: [PageNumber.TOTAL_PAGES], font: 'Arial', size: 18, color: '888888' }),
        ],
      })],
    });
  }

  // ── Portrait section: cover + inputs + chamber ────────────────────────────
  // Add a heading for Section 3 at the end of the portrait section.
  // If there are no plots, include a note here instead of making a landscape section.
  children.push(pageBreak());
  children.push(heading('3  ·  Plots', 32, 60));

  if (plots.length === 0) {
    children.push(spacer(120));
    children.push(para(
      'No plots were selected. Set the plot_* flags to true in your input file.',
      { color: '888888' }
    ));
  }

  const portraitSection = {
    properties: {
      page: {
        size:   { width: PAGE_W, height: PAGE_H },
        margin: { top: MARGIN, right: MARGIN, bottom: MARGIN, left: MARGIN },
      },
    },
    headers: { default: makeHeader() },
    footers: { default: makeFooter() },
    children,
  };

  // ── Landscape sections: one per pair of plots ─────────────────────────────
  // Group plots into pairs.  Each pair becomes its own document section with
  // landscape orientation.  Using separate sections (rather than inline page
  // breaks inside one section) is the only reliable way to switch orientation
  // mid-document in the OOXML format that docx-js targets.
  const landscapeSections = [];

  for (let p = 0; p < plots.length; p += 2) {
    const pair     = plots.slice(p, p + 2);   // 1 or 2 plots
    const pairChildren = [];

    pair.forEach((plot, localIdx) => {
      // Small spacer before the second plot on the same page
      if (localIdx === 1) pairChildren.push(spacer(200));

      pairChildren.push(heading(plot.title, 24, localIdx === 0 ? 60 : 60));

      const imgData = fs.readFileSync(plot.path);
      pairChildren.push(new Paragraph({
        spacing: { before: 40, after: 40 },
        children: [new ImageRun({
          type: 'png',
          data: imgData,
          // On a landscape page the content width is ~9 in — same as the PNG
          // dimensions, so the image fills the full width with no cropping.
          transformation: { width: IMG_W_PX, height: IMG_H_PX },
          altText: { title: plot.title, description: plot.title, name: plot.title },
        })],
      }));
    });

    landscapeSections.push({
      properties: {
        page: {
          // Pass portrait dimensions; docx-js + LANDSCAPE orientation swaps
          // them in the XML so the document renders in landscape correctly.
          size: {
            width:       PAGE_W,
            height:      PAGE_H,
            orientation: PageOrientation.LANDSCAPE,
          },
          margin: { top: MARGIN, right: MARGIN, bottom: MARGIN, left: MARGIN },
        },
      },
      headers: { default: makeHeader() },
      footers: { default: makeFooter() },
      children: pairChildren,
    });
  }

  // ── Assemble & write ──────────────────────────────────────────────────────
  const doc = new Document({
    styles: {
      default: {
        document: { run: { font: 'Arial', size: 20 } },
      },
    },
    // Portrait section first, then one landscape section per plot pair.
    sections: [portraitSection, ...landscapeSections],
  });

  const outPath = path.join(out_dir, report_name + '.docx');
  const buffer  = await Packer.toBuffer(doc);
  fs.writeFileSync(outPath, buffer);
  console.log('Report written to:', outPath);
}

main().catch(err => { console.error(err); process.exit(1); });
