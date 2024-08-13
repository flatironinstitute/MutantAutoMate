// @ts-check
// @ts-ignore
import { render, h } from "preact";
// @ts-ignore
import { useRef, useEffect } from "preact/hooks";
// @ts-ignore
import { signal, computed, useSignal, useSignalEffect } from "@preact/signals";
// @ts-ignore
import { html } from "htm/preact";
// @ts-ignore
import $3Dmol from "3dmol";

console.log($3Dmol);

const classes = {
  h2: "text-2xl",
  button:
    "bg-white text-black px-4 py-2 rounded outline outline-2 disabled:bg-gray-200 disabled:text-gray-400",
  input: "block outline outline-2 outline-gray-400 p-1",
  anchor: "text-blue-700 underline text-left",
};

const gene_name_signal = signal("NLGN1");
const residue1_signal = signal("D");
const position_signal = signal(140);
const residue2_signal = signal("Y");
const is_running_signal = signal(false);
const events_signal = signal([]);
const pdb_data_signal = signal({});

const log_signal = computed(() =>
  events_signal.value
    .map((e) => e.message)
    .filter((d) => typeof d === "string")
    .filter((d) => d !== "")
    .filter((d) => d?.length > 0)
    .join("\n")
);
const pdb_viewer_signal = signal(null);
const grantham_score_signal = computed(() => {
  return events_signal.value.find((e) => e.type === "grantham_score");
});
const charge_statement_signal = computed(() => {
  return events_signal.value.find((e) => e.type === "charge_statement");
});
const all_isoforms_signal = computed(() => {
  return events_signal.value.find((e) => e.all_isoforms)?.all_isoforms;
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
      className="mt-10 max-w-[600px] ms-auto me-auto bg-white rounded rounded-2xl p-10"
    >
      <${Inputs} />
      <${Separator} />
      <h2 className=${classes.h2}>Charge Statement</h2>
      <${Spacer} />
      <div>${charge_statement_signal?.value?.charge_statement}</div>
      <${Separator} />
      <h2 className=${classes.h2}>Grantham Score</h2>
      <${Spacer} />
      <div>${grantham_score_signal?.value?.grantham_statement}</div>
      <${Separator} />
      <h2 className=${classes.h2}>All Isoforms</h2>
      <${Spacer} />
      <div className="grid grid-cols-3">
        ${all_isoforms_signal?.value?.map((isoform) => {
          const url = `https://rest.uniprot.org/uniprotkb/${isoform}`;
          const found = sequences_signal?.value?.find(
            (d) => d.isoform === isoform
          );
          const icon = found ? "✅" : "";
          return html`<div className="flex gap-x-2">
            <span className="w-[2ch]">${icon}</span>
            <${Anchor} href=${url}>${isoform}</${Anchor}>
          </div>`;
        })}
      </div>
      <${Separator} />
      <h2 className=${classes.h2}>Filtered Isoforms</h2>
      <${Spacer} />
      <div className="grid grid-cols-4">
        ${filtered_isoforms_signal?.value?.map((isoform) => {
          const json_url = `https://rest.uniprot.org/uniprotkb/${isoform}.json`;
          const text_url = `https://rest.uniprot.org/uniprotkb/${isoform}.txt`;
          const fasta_url = `https://www.uniprot.org/uniprot/${isoform}.fasta`;
          return html`<div>${isoform}</div>
            <div>
              <${Anchor} href=${json_url}>JSON</${Anchor}>
            </div>
            <div>
              <${Anchor} href=${text_url}>Text</${Anchor}>
            </div>
            <div>
              <${Anchor} href=${fasta_url}>FASTA</${Anchor}>
            </div>`;
        })}
      </div>
      <${Separator} />
      <h2 className=${classes.h2}>PDB IDs</h2>
      <${Spacer} />
      <div>
        ${Object.entries(pdb_ids_signal?.value ?? {}).map(
          ([isoform, pdb_ids]) => {
            const isoform_pdb_links = (() => {
              if (!pdb_ids || pdb_ids.length === 0) {
                return html`<div className="text-gray-500 ml-[3ch]">
                  No PDB IDs found
                </div>`;
              } else {
                return pdb_ids.map((id) => {
                  const url = `https://www.rcsb.org/structure/${id}`;
                  return html`<div className="grid grid-cols-[repeat(3,6ch)] ml-[3ch]">
                  <div>${id}</div>
                  <${Anchor} href=${url}>RCSB</${Anchor}>
                  <button 
                  className=${classes.anchor}
                   onClick=${() => {
                     pdb_viewer_signal.value = id;
                   }}>
                   View
                  </button>
                </div>`;
                });
              }
            })();
            return html`<div>
              <div>${isoform}</div>
              ${isoform_pdb_links}
            </div>`;
          }
        )}
      </div>
      <${Separator} />
      <${PDBViewer} />
    </div>
  `;
}

function Inputs() {
  return html`
    <h2 className=${classes.h2}>MutantAutoMate</h2>
    <${Spacer} />
    <button
      className="block mb-2"
      onClick=${() => {
        gene_name_signal.value = "NLGN1";
        residue1_signal.value = "D";
        position_signal.value = 140;
        residue2_signal.value = "Y";
      }}
    >
      Example 1
    </button>
    <button
      className="block mb-2"
      onClick=${() => {
        gene_name_signal.value = "SHANK3";
        residue1_signal.value = "D";
        position_signal.value = 26;
        residue2_signal.value = "C";
      }}
    >
      Example 2
    </button>
    <${Spacer} />
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
          onInput=${(e) => (position_signal.value = e.target.value)}
        />
      </label>
      <label>
        <div>Resiude 2:</div>
        <input
          type="text"
          className=${classes.input}
          value=${residue2_signal}
          onInput=${(e) => (residue2_signal.value = e.target.value)}
        />
      </label>
    </div>

    <${Spacer} />
    <button
      className=${classes.button}
      disabled=${is_running_signal}
      onClick=${() => {
        events_signal.value = [];
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
      }}
    >
      Start
    </button>
    <${Spacer} />
    <h2>Log</h2>
    <${LogViewer} />
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

function PDBViewer() {
  const viewerDivRef = useRef(null);
  const viewersRef = useRef(null);
  // const viewer1Ref = useRef(null);
  useEffect(() => {
    if (!viewerDivRef.current) {
      return;
    }
    const viewers = $3Dmol.createViewerGrid(
      viewerDivRef.current,
      {
        rows: 2,
        cols: 1,
        control_all: true,
      }
      // { backgroundColor: "lightgrey" }
    );
    viewersRef.current = viewers;
    // console.log(viewers);
    // viewerRef.current = $3Dmol.createViewer("viewer", {
    //   defaultcolors: $3Dmol.rasmolElementColors,
    //   lowerZoomLimit: 1,
    //   upperZoomLimit: 500,
    // });
    // mutatedViewerRef.current = $3Dmol.createViewer("viewer2", {
    //   defaultcolors: $3Dmol.rasmolElementColors,
    //   lowerZoomLimit: 1,
    //   upperZoomLimit: 500,
    // });
    // viewerRef.current.setBackgroundColor("grey");
  }, []);

  useEffect(() => {
    if (!viewersRef.current || !pdb_viewer_signal.value) {
      return;
    }
    console.log(`PDB ID: ${pdb_viewer_signal.value}`);
    const viewers = viewersRef.current;
    const viewer1 = viewers[0][0];
    const viewer2 = viewers[1][0];
    viewer1.clear();
    viewer2.clear();
    (async () => {
      const pdb_data = await fetch(
        `https://files.rcsb.org/download/${pdb_viewer_signal.value}.pdb`
      ).then((res) => res.text());

      viewer1.addModel(pdb_data, "pdb");
      viewer1.setStyle({}, { cartoon: { color: "gray" } });
      viewer1.zoomTo();
      viewer1.render();

      const mutated_pdb_data = await fetch(`/mutate`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          pdb_string: pdb_data,
          residue1: residue1_signal.value,
          position: position_signal.value,
          residue2: residue2_signal.value,
        }),
      }).then((res) => res.text());

      if (mutated_pdb_data) {
        viewer2.addModel(mutated_pdb_data, "pdb");
        viewer2.setStyle({}, { cartoon: { color: "gray" } });
        viewer2.zoomTo();
        viewer2.render();
      }
    })();

    // fetch(`https://files.rcsb.org/download/${pdb_viewer_signal.value}.pdb`)
    //   .then((res) => res.text())
    //   .then((pdb_data) => {
    //     // console.log("PDB", pdb_data);

    //     viewer1.addModel(pdb_data, "pdb");
    //     viewer1.setStyle({}, { cartoon: { color: "gray" } });
    //     viewer1.zoomTo();
    //     viewer1.render();

    //     viewer2.addModel(pdb_data, "pdb");
    //     viewer2.setStyle({}, { cartoon: { color: "gray" } });
    //     viewer2.zoomTo();
    //     viewer2.render();
    //   });
  }, [pdb_viewer_signal.value]);

  return html`
    <h2 className=${classes.h2}>PDB Viewer</h2>
    <${Spacer} />
    <div>
      PDB ID:
      <span className="ml-2">${pdb_viewer_signal.value ?? `None`}</span>
    </div>
    <${Spacer} />
    <${PDBViewerInput} />
    <${Spacer} />
    <button
      className=${classes.button}
      onClick=${() => {
        const position = position_signal.value;
        const residue1 = residue1_signal.value;
        const rows = viewersRef.current;
        // const viewer = viewer1Ref.current;
        if (!rows) {
          return;
        }
        // viewer.removeAllLabels();
        // viewer.addResLabels({
        //   atom: "CA",
        //   resi: position,
        //   chain: "A",
        //   text: `${residue1}${position}${residue2}`,
        //   color: "black",
        //   fontSize: 12,
        //   backgroundColor: "white",
        // });
        // const selection = {
        //   resn: residue1,
        //   resi: position,
        // };
        // viewer.setStyle(selection, { stick: { colorscheme: "greenCarbon" } });
        // viewer.setStyle(selection, { cartoon: { color: "red" } });
        // viewer.render();
        // viewer.zoomTo(selection, 1e3);

        for (const row of rows) {
          const viewer = row[0];
          // Style for residue in red
          viewer.setStyle({ resi: position }, { cartoon: { color: "red" } });

          // Label for residue
          viewer.addLabel(
            position,
            {
              // position: "center",
              backgroundOpacity: 0.8,
              fontSize: 12,
              color: "red",
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
          viewer.zoomTo({ resi: residue1 }, 500);
        }
      }}
    >
      Zoom to: ${residue1_signal}${position_signal}
    </button>
    <${Spacer} />
    <div
      ref=${viewerDivRef}
      className="w-full h-[1000px] relative outline outline-black"
    ></div>
  `;
}

function PDBViewerInput() {
  const pdb_internal = useSignal("5OJ6");
  useSignalEffect(() => {
    if (
      pdb_viewer_signal.value &&
      pdb_internal.value !== pdb_viewer_signal.value
    ) {
      // pdb_internal.value = pdb_viewer_signal.value;
    }
  });
  return html`
    <div className="flex gap-x-4">
      <input
        type="text"
        className=${classes.input}
        id="pdb_id"
        name="pdb_id"
        value=${pdb_internal}
        onInput=${(e) => (pdb_internal.value = e.target.value)}
      />
      <button
        className=${classes.button}
        onClick=${() => {
          pdb_viewer_signal.value = pdb_internal.value;
        }}
      >
        Set
      </button>
    </div>
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
