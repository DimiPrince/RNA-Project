class EulerianPath {
    constructor() {
        this.graph = new Map();
        this.path = [];
    }

    addEdge(src, dest) {
        if (!this.graph.has(src)) {
            this.graph.set(src, []);
        }
        this.graph.get(src).push(dest);
    }

    findEulerianPath() {
        const inDegree = new Map();
        const outDegree = new Map();

        for (const [src, neighbors] of this.graph.entries()) {
            outDegree.set(src, neighbors.length);
            for (const neighbor of neighbors) {
                if (!inDegree.has(neighbor)) {
                    inDegree.set(neighbor, 0);
                }
                inDegree.set(neighbor, inDegree.get(neighbor) + 1);
            }
        }

        let startNode = null;
        let endNode = null;
        for (const node of new Set([...inDegree.keys(), ...outDegree.keys()])) {
            if ((inDegree.get(node) || 0) < (outDegree.get(node) || 0)) {
                if (!startNode) {
                    startNode = node;
                } else {
                    return null; // No Eulerian path
                }
            } else if ((inDegree.get(node) || 0) > (outDegree.get(node) || 0)) {
                if (!endNode) {
                    endNode = node;
                } else {
                    return null; // No Eulerian path
                }
            }
        }

        if (!startNode || !endNode) {
            startNode = this.graph.keys().next().value;
        }

        this.dfs(startNode);

        if (this.path.length !== [...outDegree.values()].reduce((acc, val) => acc + val, 0) + 1) {
            return null; // No Eulerian path
        }

        return this.path;
    }

    dfs(node) {
        while (this.graph.has(node) && this.graph.get(node).length > 0) {
            const nextNode = this.graph.get(node).shift();
            this.dfs(nextNode);
        }
        this.path.push(node);
    }
}

function reconstruct() {
    var gEnzymeInput = document.getElementById("gEnzymeInput").value;
    var ucEnzymeInput = document.getElementById("ucEnzymeInput").value;

    if (!gEnzymeInput || !ucEnzymeInput) {
        alert("Please enter both G-enzyme and U.C-enzyme fragments.");
        return;
    }

    var gEnzymeFragments = gEnzymeInput.split(", ").map(fragment => fragment.split(""));
    var ucEnzymeFragments = ucEnzymeInput.split(", ").map(fragment => fragment.split(""));

    var reconstructedRNA = reconstructRNA(gEnzymeFragments, ucEnzymeFragments);

    var outputDiv = document.getElementById("output");
    if (reconstructedRNA !== "No Eulerian path exists.") {
        outputDiv.innerHTML = "Reconstructed RNA sequence: " + reconstructedRNA;
    } else {
        outputDiv.innerHTML = "Error: No Eulerian path exists. Please check your input.";
    }
}

function reconstructRNA(G_enzyme_fragments, UC_enzyme_fragments) {
    const eulerianPathFinder = new EulerianPath();

    // Constructing the multigraph
    G_enzyme_fragments.forEach(fragment => {
        for (let i = 0; i < fragment.length - 1; i++) {
            eulerianPathFinder.addEdge(fragment[i], fragment[i + 1]);
        }
    });

    UC_enzyme_fragments.forEach(fragment => {
        for (let i = 0; i < fragment.length - 1; i++) {
            eulerianPathFinder.addEdge(fragment[i], fragment[i + 1]);
        }
    });

    // Finding Eulerian path and reconstructing the sequence
    const eulerianPath = eulerianPathFinder.findEulerianPath();
    if (eulerianPath) {
        const rnaSequence = eulerianPath.reverse().join('');
        return rnaSequence;
    } else {
        return "No Eulerian path exists.";
    }
}
