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
    var ucEnzymeResult = applyEnzyme(gEnzymeFragments, 'U.C');
    var gEnzymeResult = applyEnzyme(ucEnzymeFragments, 'G');
    var singleFragments = getSingleFragments(gEnzymeFragments, ucEnzymeFragments);
    var interiorExtendedBases = getInteriorExtendedBases(gEnzymeFragments, ucEnzymeFragments);
    var nonSingleFragments = getNonSingleFragments(gEnzymeFragments, ucEnzymeFragments);

    var outputDiv = document.getElementById("output");
    outputDiv.innerHTML = "<strong>Reconstructed RNA sequence:</strong> " + reconstructedRNA + "<br><br>";
    outputDiv.innerHTML += "<strong>Result of applying U.C-enzyme to G-enzyme input:</strong> " + ucEnzymeResult + "<br><br>";
    outputDiv.innerHTML += "<strong>Result of applying G-enzyme to U.C-enzyme input:</strong> " + gEnzymeResult + "<br><br>";
    outputDiv.innerHTML += "<strong>Interior extended bases from applying enzymes above:</strong> " + interiorExtendedBases + "<br><br>";
    outputDiv.innerHTML += "<strong>All non-single fragments from applying enzymes above:</strong> " + nonSingleFragments + "<br><br>";
}

function applyEnzyme(enzymeFragments, enzymeType) {
    // Function to apply the specified enzyme to the given fragments
    var enzymeResult = "";
    enzymeFragments.forEach(fragment => {
        const formattedFragment = fragment.join(""); // Joining the array to form a single string
        enzymeResult += formattedFragment + "/" + enzymeType + ", ";
    });
    return enzymeResult.slice(0, -2); // Remove the last comma and space
}

function getSingleFragments(gEnzymeFragments, ucEnzymeFragments) {
    // Function to get single fragments from applying the enzymes
    var singleFragments = [];
    gEnzymeFragments.forEach(fragment => {
        if (fragment.length === 1) {
            singleFragments.push(fragment.join(""));
        }
    });
    ucEnzymeFragments.forEach(fragment => {
        if (fragment.length === 1) {
            singleFragments.push(fragment.join(""));
        }
    });
    return singleFragments.join(" ");
}

function getInteriorExtendedBases(gEnzymeFragments, ucEnzymeFragments) {
    // Function to get interior extended bases from applying the enzymes
    var interiorExtendedBases = new Set();
    gEnzymeFragments.forEach(fragment => {
        if (fragment.length > 1) {
            for (let i = 1; i < fragment.length - 1; i++) {
                interiorExtendedBases.add(fragment[i]);
            }
        }
    });
    ucEnzymeFragments.forEach(fragment => {
        if (fragment.length > 1) {
            for (let i = 1; i < fragment.length - 1; i++) {
                interiorExtendedBases.add(fragment[i]);
            }
        }
    });
    return Array.from(interiorExtendedBases).join(" ");
}

function getNonSingleFragments(gEnzymeFragments, ucEnzymeFragments) {
    // Function to get non-single fragments from applying the enzymes
    var nonSingleFragments = [];
    gEnzymeFragments.forEach(fragment => {
        if (fragment.length > 1) {
            nonSingleFragments.push(fragment.join(""));
        }
    });
    ucEnzymeFragments.forEach(fragment => {
        if (fragment.length > 1) {
            nonSingleFragments.push(fragment.join(""));
        }
    });
    return nonSingleFragments.join(" ");
}

function reconstructRNA(G_enzyme_fragments, UC_enzyme_fragments) {
    const eulerianPathFinder = new EulerianPath();

    // Constructing the multigraph
    G_enzyme_fragments.forEach(fragment => {
        const formattedFragment = fragment.join(""); // Joining the array to form a single string
        console.log("Adding G-enzyme fragment:", formattedFragment);
        for (let i = 0; i < formattedFragment.length - 1; i++) {
            eulerianPathFinder.addEdge(formattedFragment[i], formattedFragment[i + 1]);
            console.log("Adding edge:", formattedFragment[i], "->", formattedFragment[i + 1]);
        }
    });

    UC_enzyme_fragments.forEach(fragment => {
        const formattedFragment = fragment.join(""); // Joining the array to form a single string
        console.log("Adding UC-enzyme fragment:", formattedFragment);
        for (let i = 0; i < formattedFragment.length - 1; i++) {
            eulerianPathFinder.addEdge(formattedFragment[i], formattedFragment[i + 1]);
            console.log("Adding edge:", formattedFragment[i], "->", formattedFragment[i + 1]);
        }
    });

    // Finding Eulerian path and reconstructing the sequence
    const eulerianPath = eulerianPathFinder.findEulerianPath();
    console.log("Eulerian path:", eulerianPath);
    if (eulerianPath) {
        const rnaSequence = eulerianPath.reverse().join('');
        return rnaSequence;
    } else {
        return "No Eulerian path exists.";
    }
}
