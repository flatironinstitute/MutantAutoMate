// @ts-check
// @ts-ignore
import { render, h } from "preact";
// @ts-ignore
import { useRef, useEffect } from "preact/hooks";
import {
  signal,
  computed,
  useSignal,
  // @ts-ignore
} from "@preact/signals";
// @ts-ignore
import { html } from "htm/preact";
// @ts-ignore
import $3Dmol from "3dmol";

console.log($3Dmol);

const classes = {
  h2: "text-2xl",
  h3: "text-lg font-bold",
  button:
    "bg-white text-black px-2 py-1 rounded outline outline-1 disabled:bg-gray-200 disabled:text-gray-400",
  buttonSmall: `bg-white text-black rounded outline outline-1 whitespace-nowrap disabled:bg-gray-200 disabled:text-gray-400 disabled:opacity-50 [padding-inline:1ch]`,
  input: "block outline outline-2 outline-gray-400 p-1",
  anchor: "text-blue-700 underline text-left",
};

const gene_name_signal = signal("NLGN1");
const residue1_signal = signal("D");
const position_signal = signal(140);
const residue2_signal = signal("Y");
const is_running_signal = signal(false);
const events_signal = signal([]);
const pdb_data_raw_signal = signal(null);
const pdb_data_trimmed_signal = signal(null);
const pdb_data_mutated_signal = signal(null);
const sequence_viewer_signal = signal(null);
const mutation_chain_id_signal = signal(null);
const loading_mutated_signal = signal(false);

const reset = () => {
  is_running_signal.value = false;
  events_signal.value = [];
  pdb_data_raw_signal.value = null;
  pdb_data_trimmed_signal.value = null;
  pdb_data_mutated_signal.value = null;
  sequence_viewer_signal.value = null;
  mutation_chain_id_signal.value = null;
  loading_mutated_signal.value = false;
};

const log_signal = computed(() =>
  events_signal.value
    .map((e) => e.message)
    .filter((d) => typeof d === "string")
    .filter((d) => d !== "")
    .filter((d) => d?.length > 0)
    .join("\n")
);
const grantham_score_signal = computed(() => {
  return events_signal.value.find((e) => e.type === "grantham_score");
});
const charge_statement_signal = computed(() => {
  return events_signal.value.find((e) => e.type === "charge_statement");
});
const all_isoforms_signal = computed(() => {
  return events_signal.value.find((e) => e.all_isoforms)?.all_isoforms;
});
const matching_isoforms_signal = computed(() => {
  return events_signal.value.find((e) => e.matching_isoforms)
    ?.matching_isoforms;
});
const filtered_isoforms_signal = computed(() => {
  return events_signal.value.find((e) => e.filtered_isoforms)
    ?.filtered_isoforms;
});
const pdb_ids_signal = computed(() => {
  return events_signal.value.find((e) => e.pdb_ids)?.pdb_ids;
});
const sequences_signal = computed(() => {
  return events_signal?.value?.filter((d) => d.type === "sequence") ?? [];
});

const add_event = (event) => {
  events_signal.value = [...events_signal.value, event];
};

export function App() {
  return html`
    <div
      className="mt-10 max-w-[min(1500px,90vw)] ms-auto me-auto bg-white rounded rounded-2xl p-10 gap-4 flex flex-col md:grid md:text-[0.7rem] md:grid-cols-3"
    >
      <h2 className=${classes.h2 + ` col-span-3`}>MutantAutoMate</h2>
      <div>
        <${Inputs} />
      </div>
      <div>
        <${Examples} />
      </div>
      <div className="space-y-4">
        <div>
          <h2 className=${classes.h2}>Charge Statement</h2>
          <${Spacer} />
          <div>${charge_statement_signal?.value?.charge_statement}</div>
        </div>
        <div>
          <h2 className=${classes.h2}>Grantham Score</h2>
          <${Spacer} />
          <div>${grantham_score_signal?.value?.grantham_statement}</div>
        </div>
      </div>
      <div>
        <h2 className=${classes.h2}>All Isoforms</h2>
        <${Spacer} />
        <${IsoformsTable} />
      </div>
      <div className="col-span-2">
        <h2 className=${classes.h2}>Filtered Isoforms</h2>
        <${Spacer} />
        <${IsoformCards} />
      </div>
      <div className="col-span-3">
        <h2 className=${classes.h2}>Sequence</h2>
        <${Spacer} />
        <${SequenceViewer} />
      </div>
      <div className="col-span-3">
        <${PDBViewer} />
      </div>
      <div>
        <h2 className=${classes.h2}>PDB Data</h2>
        <${Spacer} />
        <textarea
          readonly
          rows=${20}
          className="p-2 outline outline-black w-full font-mono text-xs"
          children=${pdb_data_raw_signal}
          placeholder="Raw PDB Data"
        >
        </textarea>
        <${Spacer} />
        <h2 className=${classes.h2}>Trimmed PDB</h2>
        <${Spacer} />
        <textarea
          readonly
          rows=${20}
          className="p-2 outline outline-black w-full font-mono text-xs"
          children=${pdb_data_trimmed_signal}
          placeholder="Trimmed PDB Data"
        >
        </textarea>
        <${Spacer} />
        <h2 className=${classes.h2}>Mutated PDB</h2>
        <${Spacer} />
        <${MutateButton} />
        <${Spacer} />
        <textarea
          readonly
          rows=${20}
          className="p-2 outline outline-black w-full font-mono text-xs"
          children=${pdb_data_mutated_signal}
          placeholder="Mutated PDB Data"
        >
        </textarea>
      </div>
      <div>
        <h2 className=${classes.h2}>Log</h2>
        <${Spacer} />
        <${LogViewer} />
      </div>
    </div>
  `;
}

function SequenceViewer() {
  const value = sequence_viewer_signal.value ?? "";
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
      const bg = index + 1 === position_signal.value ? "bg-yellow-200" : "";
      return html`<span className=${`relative block ${bg}`}>
        ${marker}${letter}
      </span>`;
    });
  return html`<div
    className="flex flex-wrap gap-y-[1.2rem] font-mono leading-none"
  >
    ${letters}
  </div>`;
}

async function getMutated() {
  loading_mutated_signal.value = true;
  const mutated_pdb_data = await fetch(`/mutate2`, {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({
      pdb_string: pdb_data_trimmed_signal.value,
      chain_id: mutation_chain_id_signal.value,
      position: position_signal.value,
      to_residue: residue2_signal.value,
    }),
  }).then((res) => res.text());
  pdb_data_mutated_signal.value = mutated_pdb_data;
  loading_mutated_signal.value = false;
}

function MutateButton() {
  const is_mutating_signal = useSignal(false);
  return html`
    <button
      className=${classes.button}
      disabled=${!pdb_data_trimmed_signal.value}
      onClick=${async () => {
        is_mutating_signal.value = true;
        await getMutated();
        is_mutating_signal.value = false;
      }}
    >
      Mutate PDB
    </button>
    <span className="ml-4">
      ${is_mutating_signal.value ? "Mutating..." : ""}
    </span>
    <span className="ml-4">
      ${mutation_chain_id_signal.value
        ? `Chain: ${mutation_chain_id_signal.value}`
        : ""}
    </span>
  `;
}

function IsoformCards() {
  const position = position_signal.value;
  const filtered_isoforms = filtered_isoforms_signal.value ?? [];
  const cards = filtered_isoforms.map((isoform) => {
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
    const pdb_ids = pdb_ids_signal.value?.[isoform] ?? [];
    const fetch_sequence = async () => {
      const response = await fetch(
        `https://rest.uniprot.org/uniprotkb/${isoform}.fasta`
      ).then((res) => res.text());
      sequence_viewer_signal.value = response.split("\n").slice(1).join("\n");
    };
    let pdb_rows = html`<div className="text-gray-500 ml-[3ch]">
      No PDB IDs found
    </div>`;
    if (pdb_ids.length > 0) {
      pdb_rows = pdb_ids.map(([pdb_id, chains_text]) => {
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
        const fetch_pdb = async () => {
          pdb_data_raw_signal.value = null;
          pdb_data_trimmed_signal.value = null;
          pdb_data_mutated_signal.value = null;
          // Set the mutation chain id to the first chain
          mutation_chain_id_signal.value = chains[0];
          await fetch_sequence();
          const pdb_string_raw = await fetch(
            `https://files.rcsb.org/download/${pdb_id}.pdb`
          ).then((res) => res.text());
          pdb_data_raw_signal.value = pdb_string_raw;
          const trimmed = await fetch(`/trim_pdb`, {
            method: "POST",
            headers: {
              "Content-Type": "application/json",
            },
            body: JSON.stringify({ pdb_data: pdb_string_raw, chains }),
          }).then((res) => res.text());
          pdb_data_trimmed_signal.value = trimmed;
          await getMutated();
        };
        const load_button = html`<button
          disabled=${!is_in_pdb}
          className=${classes.buttonSmall}
          onClick=${fetch_pdb}
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
          </div>
        `;
      });
    }
    const load_from_alphafold = async () => {
      await fetch_sequence();
      const url = `https://alphafold.ebi.ac.uk/api/prediction/${isoform}`;
      const response = await fetch(url).then((res) => res.json());
      const [first] = response;
      const pdb_url = first?.pdbUrl;
      const pdb_string_raw = await fetch(pdb_url).then((res) => res.text());
      pdb_data_raw_signal.value = pdb_string_raw;
      pdb_data_trimmed_signal.value = pdb_string_raw;
      await getMutated();
    };
    return html`
      <div className="space-y-2">
        <h3 className=${classes.h3}>${isoform}</h3>
        ${anchors}
        <div>PDB IDs:</div>
        <div className="space-y-2">${pdb_rows}</div>
        <button className=${classes.buttonSmall} onClick=${load_from_alphafold}>
          Load PDB from Alphafold
        </button>
      </div>
    `;
  });
  return html`<div className="space-y-4">${cards}</div>`;
}

function IsoformsTable() {
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

  return html`
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
          <th>Got Sequence</th>
          <th>Residue Match</th>
          <th>Gene Name Match</th>
        </tr>
      </thead>
      <tbody>
        ${rows}
      </tbody>
    </table>
  `;
}

function Inputs() {
  const start_processing = () => {
    reset();
    is_running_signal.value = true;
    const url = new URL("/process", window.location.origin);
    url.searchParams.append("gene_name", gene_name_signal);
    url.searchParams.append("residue1", residue1_signal);
    url.searchParams.append("position", position_signal);
    url.searchParams.append("residue2", residue2_signal);

    const eventSource = new EventSource(url);

    eventSource.onmessage = (event) => {
      console.log(event.data);
      let parsed;
      try {
        parsed = JSON.parse(event.data);
      } catch (e) {
        console.error(e);
        return;
      }
      add_event(parsed);
      if (parsed.type === "done") {
        is_running_signal.value = false;
        eventSource.close();
      }
    };

    eventSource.onerror = (error) => {
      console.error("EventSource failed:", error);
      eventSource.close();
    };
  };
  return html`
    <div className="flex flex-col gap-y-4">
      <label>
        <div>Gene Name:</div>
        <input
          type="text"
          className=${classes.input}
          value=${gene_name_signal}
          onInput=${(e) => (gene_name_signal.value = e.target.value)}
        />
      </label>
      <label>
        <div>Residue 1:</div>
        <input
          type="text"
          className=${classes.input}
          value=${residue1_signal}
          onInput=${(e) => (residue1_signal.value = e.target.value)}
        />
      </label>
      <label>
        <div>Position:</div>
        <input
          type="text"
          className=${classes.input}
          value=${position_signal}
          onInput=${(e) => (position_signal.value = +e.target.value)}
        />
      </label>
      <label>
        <div>Residue 2:</div>
        <input
          type="text"
          className=${classes.input}
          value=${residue2_signal}
          onInput=${(e) => (residue2_signal.value = e.target.value)}
        />
      </label>
      <button
        className=${classes.button + " w-[20ch]"}
        disabled=${is_running_signal}
        onClick=${start_processing}
      >
        Start
      </button>
    </div>
  `;
}

function Examples() {
  const others = `SLC6A1 S295, FRMPD4 E471, NRXN1 T324, NRXN1 V1214, ACTB V298`
    .split(`,`)
    .map((d) => d.trim())
    .filter((d) => d.length > 0)
    .map((string) => {
      const [gene_name, residue] = string.split(" ");
      const residue1 = residue[0];
      const position = +residue.slice(1);
      const residue2 = `C`;
      return [gene_name, residue1, position, residue2];
    });
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
        }}
      >
        ${gene_name} ${residue1}${position}
      </button>
    `;
  });
  return html`
    <h4 className="mb-2">Examples</h4>
    <div className="flex flex-wrap gap-4">${buttons}</div>
  `;
}

function LogViewer() {
  const textAreaRef = useRef(null);
  useEffect(() => {
    textAreaRef.current.scrollTop = textAreaRef.current.scrollHeight;
  }, [log_signal.value]);
  return html`<textarea
    ref=${textAreaRef}
    readonly
    rows=${20}
    className="p-2 outline outline-black w-full font-mono text-xs"
    children=${log_signal}
    placeholder="Log"
  >
  </textarea>`;
}

const usePDBViewer = (pdbSignal, viewer) => {
  useEffect(() => {
    if (!viewer || !pdbSignal.value) return;
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

function PDBViewer() {
  const viewerDivRef = useRef(null);
  const viewersGridRef = useRef(null);
  useEffect(() => {
    if (!viewerDivRef.current) {
      return;
    }
    viewersGridRef.current = $3Dmol.createViewerGrid(viewerDivRef.current, {
      rows: 1,
      cols: 2,
      control_all: true,
    });
  }, []);

  usePDBViewer(pdb_data_trimmed_signal, viewersGridRef.current?.[0][0]);
  usePDBViewer(pdb_data_mutated_signal, viewersGridRef.current?.[0][1]);

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
          position,
          {
            position: "center",
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
    <h2 className=${classes.h2}>PDB Viewer</h2>
    <${Spacer} />
    <button
      className=${classes.button}
      onClick=${zoom_to}
      disabled=${!pdb_data_trimmed_signal.value}
    >
      Zoom to Position: ${position_signal.value}
    </button>
    <div className="inline-block ml-4 text-lg">${loading_text}</div>
    <${Spacer} />
    <div
      ref=${viewerDivRef}
      className="w-full h-[400px] relative outline outline-black"
    ></div>
  `;
}

function Separator() {
  return html`<${Spacer} size=${10} />
    <hr className="h-[2px] bg-black border-none" />
    <${Spacer} />`;
}

function Spacer({ size = 5 }) {
  return html`<div className=${`h-${size}`}></div>`;
}

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

render(h(App), document.body);
window.scrollTo(0, 0);
