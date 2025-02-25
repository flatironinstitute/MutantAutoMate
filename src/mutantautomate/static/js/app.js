// @ts-check

/**
 * NOTE: This project **does not** use a build step for JavaScript.
 * - All of these libraries are fetched live in the browser.
 * - The `package.json` in the root is **only** used for local Typescript checks.
 * - You **do not** need to run `npm install` to deploy this page.
 * - We are using an [importmap](https://developer.mozilla.org/en-US/docs/Web/HTML/Element/script/type/importmap).
 * - The importmap is defined in the HTML file named `index2.html`.
 * - All packages are fetched from a free JS CDN, [esm.sh](https://esm.sh/).
 */
import { Fragment, render } from "preact";
import { useRef, useEffect } from "preact/hooks";
import { signal, computed, batch, effect } from "@preact/signals";
import { html } from "htm/preact";
import $3Dmol from "3dmol";
import { z } from "zod";

/** Edit this to edit the intro text. */
const intro_text = `MutantAutoMate is a tool that lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed ultricies at tortor ut facilisis. Ut commodo nibh quis nisl porttitor mattis vitae vitae augue. Nam purus mauris, accumsan sit amet vulputate a, placerat vel orci. Sed sagittis eros vel erat ullamcorper, ac faucibus est maximus. Nullam a justo non ipsum porta scelerisque. Duis mauris lacus, volutpat nec lectus ut, congue convallis arcu. Aliquam placerat massa dictum arcu pulvinar viverra et vitae justo. Morbi eu diam lorem. Vivamus sit amet vestibulum nunc, id convallis elit. Fusce augue lacus, suscipit nec mi a, commodo tincidunt metus. Integer sed dui ut nunc luctus tempus.`;
const start_message = `Pending. Click "Start" to begin processing.`;

/**
 * Top-Level Type Definitions
 * =========================================================
 */

/**
 * @template T
 * @typedef {import("@preact/signals").Signal<T>} Signal
 */

/**
 * @template T
 * @typedef {import("@preact/signals").ReadonlySignal<T>} ReadonlySignal
 */

/**
 * Schema
 * =========================================================
 */

const ShortIdSchema = z.string().min(3).max(15);
const ShortIdArraySchema = z.array(ShortIdSchema);

const MessageBaseSchema = z
  .object({
    message: z.string(),
  })
  .strict();

const MessageEventSchema = MessageBaseSchema.extend({
  type: z.literal("message"),
}).strict();

const GranthamScoreEventSchema = MessageBaseSchema.extend({
  type: z.literal("grantham_score"),
  grantham_statement: z.string(),
  grantham_score: z.number(),
}).strict();

const ChargeStatementEventSchema = MessageBaseSchema.extend({
  type: z.literal("charge_statement"),
  charge_statement: z.string(),
}).strict();

const IsoformsProgressEventSchema = MessageBaseSchema.extend({
  type: z.literal("isoforms_progress"),
  message: z.string().regex(/^Fetched \d+ isoforms.$/),
}).strict();

const AllIsoformsEventSchema = MessageBaseSchema.extend({
  type: z.literal("all_isoforms"),
  all_isoforms: ShortIdArraySchema,
}).strict();

const IsoformSequenceEventSchema = MessageBaseSchema.extend({
  type: z.literal("isoform_sequence"),
  isoform: ShortIdSchema,
}).strict();

const IsoformSequenceProgressEventSchema = MessageBaseSchema.extend({
  type: z.literal("isoform_sequence_progress"),
  message: z.string().regex(/^Fetched \d+ \/ \d+ sequences.$/),
}).strict();

const MatchingIsoformsEventSchema = MessageBaseSchema.extend({
  type: z.literal("matching_isoforms"),
  matching_isoforms: ShortIdArraySchema,
}).strict();

const GeneNameEventSchema = MessageBaseSchema.extend({
  type: z.literal("gene_name"),
  gene_name: z.string(),
}).strict();

const FilteredIsoformsEventSchema = MessageBaseSchema.extend({
  type: z.literal("filtered_isoforms"),
  filtered_isoforms: ShortIdArraySchema,
}).strict();

const PairwiseScoreSchema = z
  .object({
    isoform1: z.string(),
    isoform2: z.string(),
    score: z.number(),
  })
  .strict();

const PairwiseScoresEventSchema = MessageBaseSchema.extend({
  type: z.literal("pairwise_scores"),
  pairwise_scores: z.array(PairwiseScoreSchema),
}).strict();

const PdbIdsTupleSchema = z.tuple([
  z.string().min(3),
  z.string().min(1),
  z.string().min(1),
]);
const PdbIdsSchema = z.record(ShortIdSchema, z.array(PdbIdsTupleSchema));
const PdbIdsEventSchema = MessageEventSchema.extend({
  type: z.literal("pdb_ids"),
  pdb_ids: PdbIdsSchema,
}).strict();

const DoneEventSchema = MessageEventSchema.extend({
  type: z.literal("done"),
});

/**
 * This is the schema for any event coming from the `EventSource`
 */
const EventSchema = z
  .union([
    MessageEventSchema,
    GranthamScoreEventSchema,
    ChargeStatementEventSchema,
    IsoformsProgressEventSchema,
    AllIsoformsEventSchema,
    IsoformSequenceEventSchema,
    IsoformSequenceProgressEventSchema,
    MatchingIsoformsEventSchema,
    GeneNameEventSchema,
    FilteredIsoformsEventSchema,
    PairwiseScoresEventSchema,
    PdbIdsEventSchema,
    DoneEventSchema,
  ])
  .readonly();

/**
 * Use `zod` magic to infer different types ✨.
 *
 * @typedef {z.infer<typeof EventSchema>} Event
 */

/**
 * Tailwind CSS classes
 * =========================================================
 */

const classes = {
  h2: "text-2xl",
  h3: "text-lg font-bold",
  button:
    "bg-white text-black px-2 py-1 rounded outline outline-1 disabled:bg-gray-200 disabled:text-gray-400",
  bigButton: `text-lg bg-white text-black px-2 py-1 rounded outline outline-1 disabled:bg-gray-200 disabled:text-gray-400 w-full`,
  smallButton: `bg-white text-black rounded outline outline-1 whitespace-nowrap disabled:bg-gray-200 disabled:text-gray-400 disabled:opacity-50 [padding-inline:1ch]`,
  input: "block outline outline-2 outline-gray-400 p-1",
  bigInput: `block outline outline-2 outline-gray-400 p-1 w-full text-lg`,
  anchor: "text-blue-700 underline text-left",
};

/**
 * Signals
 * =========================================================
 */

/** @type {Signal<string>} */
const gene_name_signal = signal("NLGN1");

/** @type {Signal<string>} */
const residue1_signal = signal("D");

/** @type {Signal<number>} */
const position_signal = signal(140);

/** @type {Signal<string>} */
const residue2_signal = signal("Y");

/** @type {Signal<Event[]>} */
const events_signal = signal([
  {
    type: `message`,
    message: start_message,
  },
]);

/** @type {Signal<boolean>} */
const is_running_signal = signal(false);

/** @type {Signal<string | null>} */
const selected_isoform_signal = signal(null);

/** @type {Signal<string | null>} */
const selected_pdb_id_signal = signal(null);

/** @type {Signal<boolean>} */
const loading_mutated_signal = signal(false);

/** @type {Signal<string | null>} */
const pdb_data_trimmed_signal = signal(null);

/** @type {Signal<string | null>} */
const pdb_data_mutated_signal = signal(null);

/** @type {Signal<string | null>} */
const sequence_signal = signal(null);

/** @type {Signal<any>} */
const dssp_signal = signal(null);

effect(() => {
  const allEvents = events_signal.value;
  const latestEvent = allEvents.at(-1);
  if (!latestEvent) return;
  console.log(`Latest event:`, latestEvent);
});

const all_isoforms_signal = computed(() => {
  return events_signal.value.find((e) => e.type === "all_isoforms")
    ?.all_isoforms;
});

const sequences_signal = computed(() => {
  return (
    events_signal?.value?.filter((e) => e.type === "isoform_sequence") ?? []
  );
});

const matching_isoforms_signal = computed(() => {
  return events_signal.value.find((e) => e.type === "matching_isoforms")
    ?.matching_isoforms;
});

const filtered_isoforms_signal = computed(() => {
  return events_signal.value.find((e) => e.type === "filtered_isoforms")
    ?.filtered_isoforms;
});

const pdb_ids_signal = computed(() => {
  const pdb_ids = events_signal.value.find(
    (e) => e.type === "pdb_ids"
  )?.pdb_ids;
  return pdb_ids;
});

/** @type {ReadonlySignal<string[]>} */
const log_array_signal = computed(() =>
  events_signal.value
    .map((e) => e.message)
    .filter((d) => typeof d === "string")
    .filter((d) => d !== "")
    .filter((d) => d?.length > 0)
);

/** @type {ReadonlySignal<z.infer<typeof GranthamScoreEventSchema> | undefined>} */
const grantham_score_signal = computed(() => {
  return events_signal.value.find((e) => e.type === "grantham_score");
});

/** @type {ReadonlySignal<z.infer<typeof ChargeStatementEventSchema> | undefined>} */
const charge_statement_signal = computed(() => {
  return events_signal.value.find((e) => e.type === "charge_statement");
});

/**
 * App Code from here down
 * =========================================================
 */

/**
 * Add an event to the events array.
 * @param {Event} event
 */
function addEvent(event) {
  events_signal.value = [...events_signal.value, event];
}

/**
 * Start the `EventSource` and push events to the events signal.
 * The `EventSource`
 * @param {URL} url
 */
function startEventSource(url) {
  const eventSource = new EventSource(url);

  const onMessage = (event) => {
    try {
      /** @type {unknown} */
      let parsed;
      try {
        parsed = JSON.parse(event.data);
      } catch (error) {
        throw new Error(`Error parsing event JSON`, {
          cause: {
            error,
            event,
            eventData: event.data,
          },
        });
      }
      const parseResult = EventSchema.safeParse(parsed);
      if (!parseResult.success) {
        console.error(`Error parsing event:`, parsed);
        console.error(`Raw error:`, parseResult.error);
        console.error(`Formatted error:`, parseResult.error.format());
        // Don't throw an error here, just log it.
        // throw new Error(`Error parsing event:`, {
        //   cause: {
        //     parseResult,
        //     parsed,
        //   },
        // });
      }

      /** The event has now passed the parsing step so we assume it's valid. */
      const validEvent = /** @type {Event} */ (parsed);

      /** Add the event to the events Signal. */
      addEvent(validEvent);
      if (validEvent.type === "done") {
        is_running_signal.value = false;
        // We don't need to wait for this promise to resolve.
        autoSelectPDB();
        console.info("Done processing, closing EventSource");
        eventSource.close();
      }
    } catch (error) {
      console.error("Error processing event:", error);
      console.warn("Closing EventSource");
      eventSource.close();
    }
  };

  eventSource.addEventListener("message", onMessage);

  eventSource.addEventListener("error", (error) => {
    console.error("EventSource failed:", error);
    eventSource.removeEventListener("message", onMessage);
    eventSource.close();
  });
}

/**
 * Reset everything and start processing
 * @returns {void}
 */
function startProcessing() {
  // Reset stuff
  batch(() => {
    is_running_signal.value = true;
    events_signal.value = [];
    selected_isoform_signal.value = null;
    selected_pdb_id_signal.value = null;
    pdb_data_trimmed_signal.value = null;
    pdb_data_mutated_signal.value = null;
    sequence_signal.value = null;
    dssp_signal.value = null;
  });

  const url = new URL("/process", window.location.origin);
  url.searchParams.append("gene_name", gene_name_signal.value);
  url.searchParams.append("residue1", residue1_signal.value);
  url.searchParams.append("position", position_signal.value.toString());
  url.searchParams.append("residue2", residue2_signal.value);

  startEventSource(url);
}

/**
 * @param {string} isoform
 * @returns {Promise<void>}
 */
async function fetchSequence(isoform) {
  console.log(`Fetching sequence for ${isoform}`);
  const response = await fetch(
    `https://rest.uniprot.org/uniprotkb/${isoform}.fasta`
  ).then((res) => res.text());
  sequence_signal.value = response.split("\n").slice(1).join("\n");
}

/**
 * @param {{ pdb_string: string }} params
 * @returns {Promise<void>}
 */
async function getDSSP({ pdb_string }) {
  const dssp_data = await fetch(`/dssp`, {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({
      pdb_string,
    }),
  }).then((res) => res.json());
  dssp_signal.value = dssp_data;
}

/**
 * @param {{ pdb_string: string; chain_id: string | null; position: number; to_residue: string }} params
 * @returns {Promise<void>}
 */
async function getMutated({ pdb_string, chain_id, position, to_residue }) {
  loading_mutated_signal.value = true;
  const mutated_pdb_data = await fetch(`/mutate`, {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({
      pdb_string,
      chain_id,
      position,
      to_residue,
    }),
  }).then((res) => res.text());
  pdb_data_mutated_signal.value = mutated_pdb_data;
  loading_mutated_signal.value = false;
}

/**
 * @param {{ isoform: string; pdb_id: string; chains: string[] | null; }} params
 * @returns {Promise<void>}
 */
async function fetchPDB({ isoform, pdb_id, chains }) {
  console.log(`Fetching PDB:`, { isoform, pdb_id, chains });
  selected_isoform_signal.value = isoform;
  selected_pdb_id_signal.value = pdb_id;
  pdb_data_trimmed_signal.value = null;
  pdb_data_mutated_signal.value = null;
  dssp_signal.value = null;
  await fetchSequence(isoform);
  const pdb_string_raw = await fetch(
    `https://files.rcsb.org/download/${pdb_id}.pdb`
  ).then((res) => res.text());
  // NOTE: Just trim to the first chain
  const first_chain = chains?.[0];
  if (!first_chain) {
    throw new Error(`No chains found for ${isoform}`);
    return;
  }
  console.log(`Trimming PDB to chain: ${first_chain}`);
  const trimmed = await fetch(`/trim_pdb`, {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({ pdb_data: pdb_string_raw, chains: [first_chain] }),
  }).then((res) => res.text());
  pdb_data_trimmed_signal.value = trimmed;
  await getDSSP({ pdb_string: trimmed });
  await getMutated({
    pdb_string: trimmed,
    chain_id: first_chain,
    position: position_signal.value,
    to_residue: residue2_signal.value,
  });
}

/**
 * @param {string} isoform
 * @returns {Promise<void>}
 */
async function fetchPDBFromAlphaFold(isoform) {
  console.log(`Fetching from AlphaFold: ${isoform}`);
  selected_isoform_signal.value = isoform;
  selected_pdb_id_signal.value = `AlphaFold`;
  pdb_data_trimmed_signal.value = null;
  pdb_data_mutated_signal.value = null;
  dssp_signal.value = null;
  if (!isoform) return;
  await fetchSequence(isoform);
  const fixed_isoform = isoform.split("-")[0];
  const url = `https://alphafold.ebi.ac.uk/api/prediction/${fixed_isoform}`;
  const response = await fetch(url).then((res) => res.json());
  const [first] = response;
  const pdb_url = first?.pdbUrl;
  const pdb_string_raw = await fetch(pdb_url).then((res) => res.text());
  // We do not need to trim
  pdb_data_trimmed_signal.value = pdb_string_raw;
  await getDSSP({ pdb_string: pdb_string_raw });
  await getMutated({
    pdb_string: pdb_string_raw,
    chain_id: null,
    position: position_signal.value,
    to_residue: residue2_signal.value,
  });
}

/**
 * Try to auto-select a PDB based on the first isoform.
 *
 * @returns {Promise<void>}
 */
async function autoSelectPDB() {
  const position = position_signal.value;
  const filtered_isoforms = filtered_isoforms_signal.value ?? [];
  const all_pdb_ids = pdb_ids_signal.value ?? {};
  console.log(`Auto-Selecting PDB`, { filtered_isoforms });
  if (filtered_isoforms.length === 0) return;

  const first_isoform = filtered_isoforms[0];

  const first_isoform_pdb_ids = all_pdb_ids[first_isoform] ?? [];

  // Auto-select a PDB. The PDB must be:
  // 1. In the first isoform's list of PDBs
  // 2. The position must be within the range of the PDB
  // 3. We prefer the PDB with the highest resolution

  const valid_pdbs = [];

  for (const [pdb_id, chains_text, resolution_text] of first_isoform_pdb_ids) {
    let chains = null;
    let is_in_pdb = false;
    let position_range = null;
    const resolution = parseFloat(resolution_text);
    if (chains_text) {
      const split = chains_text.split("=");
      const chain_letters = split[0];
      position_range = split[1];
      chains = chain_letters.split("/");
      const [position_start, position_end] = position_range.split("-");
      is_in_pdb = +position >= +position_start && +position <= +position_end;
    }
    if (is_in_pdb) {
      valid_pdbs.push({ pdb_id, chains, resolution });
    }
  }

  let found_pdb = null;

  if (valid_pdbs.length === 0) {
    // Fetch from AlphaFold
    found_pdb = {
      pdb_id: `alphafold`,
      chains: null,
      resolution: NaN,
    };
  } else {
    // Get the highest resolution PDB
    let highest_resolution = Infinity;
    for (const pdb of valid_pdbs) {
      if (pdb.resolution < highest_resolution) {
        highest_resolution = pdb.resolution;
        found_pdb = pdb;
      }
    }
  }

  if (!found_pdb) {
    throw new Error(`found_pdb is undefined for some reason.`);
  }
  const { pdb_id, chains } = found_pdb;
  if (pdb_id === `alphafold`) {
    await fetchPDBFromAlphaFold(first_isoform);
  } else {
    await fetchPDB({ isoform: first_isoform, pdb_id, chains });
  }
}

/** @returns {preact.VNode} */
function Spacer({ size = 5 }) {
  return html`<div className=${`h-${size}`}></div>`;
}

/** @returns {preact.VNode} */
function Anchor({ children, href }) {
  return html`<a
    href=${href}
    target="_blank"
    rel="noreferrer"
    className=${classes.anchor}
  >
    ${children}
  </a>`;
}

/** @returns {preact.VNode} */
function Inputs() {
  return html`
    <div class="grid grid-cols-2 gap-y-4 gap-x-4">
      <label>
        <div>Gene Name</div>
        <input
          type="text"
          className=${classes.bigInput}
          value=${gene_name_signal}
          onInput=${(e) => (gene_name_signal.value = e.target.value)}
        />
      </label>
      <label>
        <div>Residue 1</div>
        <input
          type="text"
          className=${classes.bigInput}
          value=${residue1_signal}
          onInput=${(e) => (residue1_signal.value = e.target.value)}
        />
      </label>
      <label>
        <div>Position</div>
        <input
          type="text"
          className=${classes.bigInput}
          value=${position_signal}
          onInput=${(e) => (position_signal.value = +e.target.value)}
        />
      </label>
      <label>
        <div>Residue 2</div>
        <input
          type="text"
          className=${classes.bigInput}
          value=${residue2_signal}
          onInput=${(e) => (residue2_signal.value = e.target.value)}
        />
      </label>
      <div class="col-span-2">
        <button
          className=${classes.bigButton}
          disabled=${is_running_signal}
          onClick=${startProcessing}
          data-test-id="start-button"
        >
          Start
        </button>
      </div>
      <div className="col-span-2">
        <${StatusDisplay} />
      </div>
    </div>
  `;
}

/** @returns {preact.VNode} */
function StatusDisplay() {
  const latestLog = log_array_signal.value.at(-1) ?? start_message;

  console.log(`Latest Log:`, latestLog);

  return html`
    <${Fragment}>
      <div>Status:</div>
      <div>${log_array_signal.value.at(-1)}</div>
      <div>
        <${ProgressBar} string=${latestLog} />
      </div>
    </${Fragment}>
  `;
}

/**
 * @param {{ string: string }} props
 * @returns {preact.VNode}
 */
function ProgressBar({ string }) {
  let showProgress = false;
  let progressValue = 0;
  let progressMax = 100;

  const progressRegex = /(?<first>\d+)\s\/\s(?<second>\d+)/;

  const match = progressRegex.exec(string);

  if (match !== null) {
    if (match.groups) {
      const { first, second } = match.groups;
      showProgress = true;
      progressValue = +first;
      progressMax = +second;
    }
  }

  return html`
    <${Fragment}>
      <progress
        value=${progressValue}
        max=${progressMax}
        style=${{
          opacity: showProgress ? 1 : 0,
        }}
      ></progress>
      <style>
        progress {
          --border-radius: 4px;
          width: 100%;
          height: 20px;
          margin-top: 1rem;
          -webkit-appearance: none;
          appearance: none;
          border-radius: var(--border-radius);
          overflow: hidden;
        }
        progress::-webkit-progress-bar {
          background-color: lightgrey;
          border-radius: var(--border-radius);
        }
        progress::-webkit-progress-value {
          background-color: var(--ccb-green);
          border-radius: var(--border-radius);
        }
        progress::-moz-progress-bar {
          background-color: var(--ccb-green);
          border-radius: var(--border-radius);
        }
      </style>
    </${Fragment}>
  `;
}

/** @returns {preact.VNode} */
function DSSPText() {
  const selected_position = position_signal.value;
  const dssp_data = dssp_signal.value ?? [];
  const dssp_description_string = [];
  const textMap = {
    H: "Helix",
    E: "Strand",
    C: "Coil",
  };
  for (const [frame_index, frame] of dssp_data.entries()) {
    if (frame_index > 0) continue;
    for (const [position, assignment] of frame.entries()) {
      if (position === selected_position) {
        dssp_description_string.push(
          `The secondary structure at position ${position} is "${textMap[assignment]}"`
        );
      }
    }
  }
  let text = html`Pending`;
  if (dssp_description_string.length > 0) {
    text = html`${dssp_description_string.map((d) => html`<div>${d}</div>`)}`;
  }
  return html`
    <div>
      <h2 className=${classes.h2}>DSSP</h2>
      <${Spacer} />
      <div>${text}</div>
    </div>
  `;
}

/** @returns {preact.VNode} */
function TextBoxes() {
  const isoform = selected_isoform_signal.value;
  const pdb_id = selected_pdb_id_signal.value;

  const isoform_anchors = isoform
    ? html`<div><${IsoformAnchors} isoform=${isoform} /></div>`
    : null;
  const pdb_url = `https://www.rcsb.org/structure/${pdb_id}`;
  const pdb_anchor = pdb_id
    ? html`<${Anchor} href=${pdb_url}>${pdb_id}</${Anchor}>`
    : `Pending`;
  return html`
    <div className="grid grid-cols-4 gap-4 px-10">
      <div>
        <h2 className=${classes.h2}>Selected Isoform</h2>
        <div>${selected_isoform_signal.value ?? `Pending`}</div>
        ${isoform_anchors}
        <${Spacer} />
        <h2 className=${classes.h2}>Selected PDB</h2>
        <div>${pdb_anchor}</div>
      </div>
      <div>
        <h2 className=${classes.h2}>Charge Statement</h2>
        <${Spacer} />
        <div>
          ${charge_statement_signal?.value?.charge_statement ?? `Pending`}
        </div>
      </div>
      <div>
        <h2 className=${classes.h2}>Grantham Score</h2>
        <${Spacer} />
        <div>
          ${grantham_score_signal?.value?.grantham_statement ?? `Pending`}
        </div>
      </div>
      <${DSSPText} />
    </div>
  `;
}

/**
 * @param {Signal<string | null>} pdbSignal
 * @param {$3Dmol.GLViewer | undefined} viewer
 */
const usePDBViewer = (pdbSignal, viewer) => {
  useEffect(() => {
    if (!viewer) return;
    if (!pdbSignal.value) {
      viewer.clear();
      return;
    }
    const pdb_data = pdbSignal.value;
    viewer.clear();
    (async () => {
      viewer.addModel(pdb_data, "pdb");
      viewer.setStyle({}, { cartoon: { color: "gray" } });
      viewer.zoomTo();
      viewer.render();
    })();
  }, [pdbSignal.value]);
};

/**
 * @returns {preact.VNode}
 */
function PDBViewer() {
  /** @typedef {$3Dmol.GLViewer[][]} ViewerGrid */

  const viewerDivRef = useRef(/** @type {HTMLDivElement | null} */ (null));

  const viewersGridRef = useRef(/** @type {ViewerGrid | null} **/ (null));

  useEffect(() => {
    if (!viewerDivRef.current) {
      return;
    }
    const viewerGrid = $3Dmol.createViewerGrid(viewerDivRef.current, {
      rows: 1,
      cols: 2,
      control_all: true,
    });
    viewersGridRef.current = viewerGrid;
  }, []);

  /** @type {$3Dmol.GLViewer | undefined} */
  const viewer1 = viewersGridRef.current?.[0][0];
  /** @type {$3Dmol.GLViewer | undefined} */
  const viewer2 = viewersGridRef.current?.[0][1];

  usePDBViewer(pdb_data_trimmed_signal, viewer1);
  usePDBViewer(pdb_data_mutated_signal, viewer2);

  const zoom_to = () => {
    const position = +(position_signal?.value ?? 0);
    const grid = viewersGridRef.current;
    if (!grid) return;

    for (const row of grid) {
      for (const viewer of row) {
        // Style for residue in red
        viewer.setStyle(
          { resi: position },
          { stick: { color: "red" }, cartoon: { color: "green", opacity: 0.5 } }
        );
        // Label for residue
        viewer.addLabel(
          position.toString(),
          {
            backgroundOpacity: 0.8,
            fontSize: 12,
            showBackground: true,
          },
          { resi: position }
        );
        // Additional styling and configurations
        // viewer.setStyle(
        //   { hetflag: true },
        //   { stick: { colorscheme: "greenCarbon", radius: 0.25 } }
        // ); // Heteroatoms
        // viewer.setStyle({ bonds: 0 }, { sphere: { radius: 0.5 } }); // Water molecules
        // Zoom and render
        viewer.render();
        // viewer.zoomTo({ resi: position, chain: "A" }, 500);
        viewer.zoomTo({ resi: position }, 500);
      }
    }
  };

  const loading_text = loading_mutated_signal.value
    ? "Loading mutated PDB..."
    : "";

  return html`
    <div class="px-10">
      <h2 className=${classes.h2}>PDB Viewer</h2>
      <${Spacer} />
      <button
        className=${classes.button}
        onClick=${zoom_to}
        disabled=${!pdb_data_trimmed_signal.value}
      >
        Zoom to Position: ${position_signal.value}
      </button>
      <div className="inline-block ml-4 text-xl font-bold">${loading_text}</div>
      <${Spacer} />
      <div
        ref=${viewerDivRef}
        className="w-full h-[400px] relative outline outline-black"
      ></div>
    </div>
  `;
}

/**
 * @returns {preact.VNode}
 */
function SequenceViewer() {
  const value = sequence_signal.value ?? "";
  const letters = value
    .split("")
    .filter((letter) => letter.match(/[A-Z]/))
    .map((letter, index) => {
      let marker = null;
      if ((index + 1) % 50 === 0) {
        marker = html`<span
          className="text-xs text-gray-500 absolute leading-[0] -top-[7px]"
        >
          ${index + 1}
        </span>`;
      }
      const bg =
        index + 1 === position_signal.value
          ? "bg-yellow-200 outline outline-2 outline-red-500"
          : "";
      return html`<span className=${`relative block ${bg}`}>
        ${marker}${letter}
      </span>`;
    });
  const seq = html`<div
    className="flex flex-wrap gap-y-[1.2rem] font-mono leading-none"
  >
    ${letters}
  </div>`;
  return html`
    <div class="px-10">
      <h2 class=${classes.h2}>Sequence</h2>
      <${Spacer} />
      ${value === "" ? `Pending` : seq}
    </div>
  `;
}

/**
 * @param {{ isoform: string }} props
 * @returns {preact.VNode}
 */
function IsoformAnchors({ isoform }) {
  const uniprot_url = `https://www.uniprot.org/uniprotkb/${isoform}`;
  const text_url = `https://rest.uniprot.org/uniprotkb/${isoform}.txt`;
  const json_url = `https://rest.uniprot.org/uniprotkb/${isoform}.json`;
  const fasta_url = `https://rest.uniprot.org/uniprotkb/${isoform}.fasta`;
  const uniprot_anchor = html`<${Anchor} href=${uniprot_url}>UniProt</${Anchor}>`;
  const text_anchor = html`<${Anchor} href=${text_url}>Text</${Anchor}>`;
  const json_anchor = html`<${Anchor} href=${json_url}>JSON</${Anchor}>`;
  const fasta_anchor = html`<${Anchor} href=${fasta_url}>FASTA</${Anchor}>`;
  // prettier-ignore
  const anchors = html`${uniprot_anchor} | ${text_anchor} | ${json_anchor} | ${fasta_anchor}`;
  return anchors;
}

/**
 * @returns {preact.VNode}
 */
function MatchingIsoforms() {
  const position = position_signal.value;
  const filtered_isoforms = filtered_isoforms_signal.value ?? [];
  const cards = filtered_isoforms.map((isoform) => {
    const anchors = html`<${IsoformAnchors} isoform=${isoform} />`;
    const pdb_ids = pdb_ids_signal.value?.[isoform] ?? [];
    let pdb_rows = html`<div className="text-gray-500 ml-[3ch]">
      No PDB IDs found
    </div>`;
    if (pdb_ids.length > 0) {
      const rows = pdb_ids.map(([pdb_id, chains_text, resolution_text]) => {
        let chains = null;
        let is_in_pdb = false;
        let position_range = null;
        if (chains_text) {
          const split = chains_text.split("=");
          const chain_letters = split[0];
          position_range = split[1];
          chains = chain_letters.split("/");
          const [position_start, position_end] = position_range.split("-");
          is_in_pdb =
            +position >= +position_start && +position <= +position_end;
        }
        const url = `https://www.rcsb.org/structure/${pdb_id}`;
        const load_button = html`<button
          disabled=${!is_in_pdb}
          className=${classes.smallButton}
          onClick=${() => {
            fetchPDB({ isoform, pdb_id, chains });
          }}
        >
          Load PDB
        </button>`;
        const chains_list = chains ? chains.join(`, `) : `-`;
        const pdb_anchor = html`<${Anchor} href=${url}>${pdb_id}</${Anchor}>`;
        return html`
          <div className="ml-[3ch] flex gap-x-2">
            <span className="font-mono w-[5ch]">${pdb_anchor}</span>
            ${load_button}
            <span>Chains: ${chains_list}</span>
            <span>Range: ${position_range}</span>
            <span>Resolution: ${resolution_text}</span>
          </div>
        `;
      });
      pdb_rows = html`${rows}`;
    }
    return html`
      <div className="space-y-2">
        <h3 className=${classes.h3}>${isoform}</h3>
        ${anchors}
        <div>PDB IDs:</div>
        <div className="space-y-2">${pdb_rows}</div>
        <button
          className=${classes.smallButton}
          onClick=${() => {
            fetchPDBFromAlphaFold(isoform);
          }}
        >
          Load PDB from Alphafold
        </button>
      </div>
    `;
  });
  return html` <div class="px-10">
    <h2 class=${classes.h2}>Matching Isoforms</h2>
    <${Spacer} />
    <div className="space-y-4">
      ${filtered_isoforms.length === 0 ? `Pending` : cards}
    </div>
  </div>`;
}

/**
 * @returns {preact.VNode}
 */
function AllIsoforms() {
  const all_isoforms = all_isoforms_signal.value ?? [];
  const sequences = sequences_signal.value ?? [];
  const residue_matches = matching_isoforms_signal.value ?? [];
  const gene_name_matches = filtered_isoforms_signal.value ?? [];
  const rows = all_isoforms.map((isoform) => {
    // const url = `https://rest.uniprot.org/uniprotkb/${isoform}`;
    const url = `https://www.uniprot.org/uniprotkb/${isoform}`;
    const anchor = html`<${Anchor} href=${url}>${isoform}</${Anchor}>`;
    const found_sequence = sequences.find((d) => d.isoform === isoform);
    const found_residue = residue_matches.find((d) => d === isoform);
    const found_gene_name = gene_name_matches.find((d) => d === isoform);
    return html`
      <tr>
        <td>${anchor}</td>
        <td className="text-center">${found_sequence ? `✅` : `-`}</td>
        <td className="text-center">${found_residue ? `✅` : `-`}</td>
        <td className="text-center">${found_gene_name ? `✅` : `-`}</td>
      </tr>
    `;
  });

  const table = html`
    <style>
      table.grid {
        thead,
        tbody,
        tr {
          display: contents;
        }
      }
    </style>
    <table className="grid grid-cols-4">
      <thead>
        <tr className="text-xs">
          <th>Isoform</th>
          <th>Sequence?</th>
          <th>Residue?</th>
          <th>Gene Name?</th>
        </tr>
      </thead>
      <tbody>
        ${rows}
      </tbody>
    </table>
  `;

  return html`<div class="px-10">
    <h2 class=${classes.h2}>All Isoforms</h2>
    <${Spacer} />
    ${table}
  </div>`;
}

/**
 * @returns {preact.VNode}
 */
function Examples() {
  /** @typedef {[string, string, number, string]} ExampleTuple  */
  const others = `SLC6A1 S295, FRMPD4 E471, NRXN1 T324, NRXN1 V1214, ACTB V298`
    .split(`,`)
    .map((d) => d.trim())
    .filter((d) => d.length > 0)
    .map((string) => {
      const [gene_name, residue] = string.split(" ");
      const residue1 = residue[0];
      const position = +residue.slice(1);
      const residue2 = `C`;
      /** @type {ExampleTuple} */
      const tuple = [gene_name, residue1, position, residue2];
      return tuple;
    });
  /** @type {ExampleTuple[]} */
  const examples = [
    [`NLGN1`, `D`, 140, `Y`],
    [`SHANK3`, `D`, 26, `Y`],
    [`NRXN1`, `K`, 287, `E`],
    [`CSNK1G1`, `T`, 140, `C`],
    [`SCN2A`, `R`, 1635, `C`],
    ...others,
  ];
  const buttons = examples.map(([gene_name, residue1, position, residue2]) => {
    return html`
      <button
        className=${classes.button}
        onClick=${() => {
          gene_name_signal.value = gene_name;
          residue1_signal.value = residue1;
          position_signal.value = position;
          residue2_signal.value = residue2;
          startProcessing();
        }}
      >
        ${gene_name} ${residue1}${position}
      </button>
    `;
  });
  return html`
    <section class="px-10">
      <h2 class=${classes.h2}>Examples</h2>
      <${Spacer} />
      <div className="flex flex-wrap gap-4">${buttons}</div>
    </section>
  `;
}

/**
 * @returns {preact.VNode}
 */
function App() {
  return html`
    <main class="space-y-4">
      <div className="grid lg:grid-cols-2">
        <section
          class="py-10 px-10"
          style=${{ backgroundColor: `var(--ccb-green)` }}
        >
          <h1 class="text-[3rem]">MutantAutoMate</h1>
          <p>${intro_text}</p>
        </section>
        <section class="py-10 px-10">
          <h2 class="text-3xl">Inputs</h2>
          <${Inputs} />
        </section>
      </div>
      <${TextBoxes} />
      <${PDBViewer} />
      <${SequenceViewer} />
      <${Spacer} />
      <${MatchingIsoforms} />
      <${Spacer} />
      <${AllIsoforms} />
      <${Spacer} />
      <${Examples} />
      <${Spacer} />
    </main>
  `;
}

render(html`<${App} />`, document.body);
window.scrollTo(0, 0);

// Auto-start, useful when testing
// document.addEventListener("DOMContentLoaded", () => {
//   /** @type {HTMLButtonElement | null} */
//   const startButton = document.querySelector(`[data-test-id="start-button"]`);
//   if (startButton) {
//     startButton.click();
//   }
// });
