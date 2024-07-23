// @ts-check
// @ts-ignore
import { render, h } from "preact";
// @ts-ignore
import { useReducer, useState, useRef, useEffect } from "preact/hooks";
// @ts-ignore
import { signal } from "@preact/signals";
// @ts-ignore
import { html } from "htm/preact";
// @ts-ignore
import $3Dmol from "3dmol";

console.log($3Dmol);

const classes = {
  h2: "text-2xl",
  button: "bg-white text-black px-4 py-2 rounded outline outline-2",
  input: "block outline outline-2 outline-gray-400 p-1",
};

export function App() {
  return html`
    <div
      className="mt-10 max-w-[500px] ms-auto me-auto bg-white rounded rounded-2xl p-10"
    >
      <${SearchTest} />
      <${Separator} />
      <${StreamTest} />
      <${Separator} />
      <${GranthamScore} />
      <${Separator} />
      <${PDBViewer} />
    </div>
  `;
}

function SearchTest() {
  const [geneName, setGeneName] = useState("NLGN1");
  const [output, setOutput] = useState("this will be output");

  return html`
    <h2 className=${classes.h2}>Search Test</h2>
    <${Spacer} />
    <label for="aa1">Gene Name:</label>
    <input
      type="text"
      className=${classes.input}
      id="aa1"
      name="aa1"
      value=${geneName}
      onInput=${(e) => setGeneName(e.target.value)}
    />
    <${Spacer} />
    <button
      className=${classes.button}
      onClick=${() => {
        const url = new URL("/search-uniprot", window.location.origin);
        url.searchParams.append("gene_name", geneName);
        fetch(url, {
          method: "GET",
          headers: {
            "Content-Type": "application/json",
          },
        })
          .then((res) => res.json())
          .then((data) => {
            setOutput(JSON.stringify(data, null, 2));
          });
      }}
    >
      Search
    </button>
    <${Spacer} />
    <textarea
      id="score"
      name="score"
      readonly
      rows=${10}
      className="outline outline-black w-full font-mono text-xs"
      children=${output}
    >
    </textarea>
  `;
}

function StreamTest() {
  const [output, setOutput] = useState("this will be output");
  return html`
    <h2 className=${classes.h2}>Stream Test</h2>
    <${Spacer} />
    <button
      className=${classes.button}
      onClick=${() => {
        const eventSource = new EventSource("/stream-test");
      }}
    >
      Start Stream Test
    </button>
    <${Spacer} />
    <textarea
      id="score"
      name="score"
      readonly
      rows=${10}
      className="outline outline-black w-full font-mono text-xs"
      children=${output}
    >
    </textarea>
  `;
}

function GranthamScore() {
  const [aa1, setAA1] = useState("D");
  const [aa2, setAA2] = useState("Y");
  const [result, setResult] = useState("");

  return html`
    <h2 className="text-2xl">Grantham Score</h2>
    <${Spacer} />
    <div>
      <label for="aa1">Amino Acid 1:</label>
      <input
        type="text"
        className=${classes.input}
        id="aa1"
        name="aa1"
        value=${aa1}
        onInput=${(e) => setAA1(e.target.value)}
      />
    </div>
    <${Spacer} />
    <div>
      <label for="aa2">Amino Acid 2:</label>
      <input
        type="text"
        className=${classes.input}
        id="aa2"
        name="aa2"
        value=${aa2}
        onInput=${(e) => setAA2(e.target.value)}
      />
    </div>
    <${Spacer} />
    <button
      className=${classes.button}
      onClick=${(e) => {
        const url = new URL("/grantham", window.location.origin);
        url.searchParams.append("amino_acid_1", aa1);
        url.searchParams.append("amino_acid_2", aa2);
        fetch(url, {
          method: "GET",
          headers: {
            "Content-Type": "application/json",
          },
        })
          .then((res) => res.json())
          .then((data) => {
            setResult(JSON.stringify(data, null, 2));
          });
      }}
    >
      Get Grantham Score
    </button>
    <${Spacer} />
    <div>
      <textarea
        id="score"
        name="score"
        readonly
        rows=${10}
        className="outline outline-black w-full font-mono text-xs"
        children=${result}
      >
      </textarea>
    </div>
  `;
}

function PDBViewer() {
  const [pdbId, setPdbId] = useState("1BNA");
  const viewerRef = useRef(null);
  useEffect(() => {
    viewerRef.current = $3Dmol.createViewer("viewer", {
      defaultcolors: $3Dmol.rasmolElementColors,
    });
    viewerRef.current.setBackgroundColor(0xff00ff);
  }, []);
  return html`
    <h2 className=${classes.h2}>PDB Viewer</h2>
    <${Spacer} />
    <div>
      <label>PDB ID:</label>
      <input
        type="text"
        className=${classes.input}
        id="pdb_id"
        name="pdb_id"
        value=${pdbId}
      />
    </div>
    <${Spacer} />
    <button
      className=${classes.button}
      onClick=${(e) => {
        const viewer = viewerRef.current;
        viewer.clear();
        viewer.setBackgroundColor(0xffffff);
        $3Dmol.download(`pdb:${pdbId}`, viewer, {}, () => {
          viewer.setStyle({}, { cartoon: { color: "spectrum" } });
          viewer.zoomTo();
          viewer.render();
        });
        // const url = new URL(
        //   `https://files.rcsb.org/view/${pdbId}.pdb`,
        //   window.location.origin
        // );
        // viewer.current.load("rcsb://" + pdbId, "pdb", () => {
        //   viewer.current.zoomTo();
        //   viewer.current.render();
        // fetch(url)
        //   .then((response) => {
        //     if (!response.ok) {
        //       throw new Error(`HTTP error! status: ${response.status}`);
        //     }
        //     return response.text();
        //   })
        //   .then((data) => {
        //     viewer.addModel(data, "pdb"); // Load the new PDB data
        //     viewer.setStyle({}, { cartoon: { color: "spectrum" } }); // Set style
        //     viewer.zoomTo(); // Zoom to fit the model
        //     viewer.render(); // Render the model
        //   })
        //   .catch((e) => {
        //     console.error(e);
        //   });
      }}
    >
      Load
    </button>
    <${Spacer} />
    <div
      id="viewer"
      className="w-full h-[500px] relative outline outline-black"
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

render(h(App), document.body);
