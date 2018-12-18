'use strict';

const HORIZONTAL_SCALE = 100;

class Graph {
    constructor(nodes, edges) {
        this.nodes = {};
        for (let node of nodes) {
            if (!node.dummy) {
                this.nodes[node.id] = node;
            }
        }
        this.edges = {};
        for (let node of nodes) {
            this.edges[node.id] = [];
        }
        
        for (let edge of edges) {
            if (edge.source.dummy || edge.target.dummy) {
                continue;
            }

            const v1 = edge.source.id;
            const v2 = edge.target.id;
            if (edge.left) {
                this.edges[v2].push(v1);
            } else if (edge.right) {
                this.edges[v1].push(v2);
            }
        }
    }
}

function sugiyama(nodes, edges, withDummy) {
    const graph = new Graph(nodes, edges);

    if (Object.keys(graph.nodes) == 0) {
        return [nodes, edges];
    }

    const layers = longestPathLayering(graph);
    const maxLayer = Math.max(...Object.values(layers));

    addDummyVertices(graph, layers);

    let nodesByLayer = new Array(maxLayer + 1).fill().map(() => []);
    for (let vertex in layers) {
        nodesByLayer[layers[vertex]].push(parseInt(vertex));
    }

    nodesByLayer = edgeCrossingMinimization(graph, nodesByLayer);

    computeHorizontalPlacement(graph, nodesByLayer);

    for (let layer = 0; layer < nodesByLayer.length; ++layer) {
        for (let i = 0; i < nodesByLayer[layer].length; ++i) {
            graph.nodes[nodesByLayer[layer][i]].y = 50 + 100 * layer;
            graph.nodes[nodesByLayer[layer][i]].x = width / 2 + graph.nodes[nodesByLayer[layer][i]].x;
        }
    }

    return convertToEditorFormat(graph, withDummy);
}

function dfs(graph, visited, topsort, vertex) {
    if (visited[vertex] == 2) {
        return;
    } else if (visited[vertex] == 1) {
        throw new Error('Graph should be acyclic');
    }
    visited[vertex] = 1;
    for (let neighbor of graph.edges[vertex]) {
        dfs(graph, visited, topsort, neighbor);
    }
    visited[vertex] = 2;
    topsort.push(vertex);
}

function longestPathLayering(graph) {
    const visited = {};
    const topsort = [];

    for (let vertex of Object.keys(graph.nodes)) {
        visited[vertex] = 0;
    }

    for (let vertex of Object.keys(graph.nodes)) {
        if (visited[vertex] == 0) {
            dfs(graph, visited, topsort, vertex);
        }
    }

    const layers = {};

    for (let node of Object.values(graph.nodes)) {
        layers[node.id] = 0;
    }

    for (let vertex of topsort.reverse()) {
        for (let neighbor of graph.edges[vertex]) {
            layers[neighbor] = Math.max(layers[neighbor], layers[vertex] + 1);
        }
    }

    for (let vertex of topsort.reverse()) {
        if (graph.edges[vertex].length > 0) {
            layers[vertex] = layers[graph.edges[vertex][0]] - 1;
        }
        for (let neighbor of graph.edges[vertex]) {
            layers[vertex] = Math.min(layers[vertex], layers[neighbor] - 1);
        }
    }

    return layers;
}

function addDummyVertices(graph, layers) {
    let nextAvailableDummyId = -1;
    for (let node of Object.values(graph.nodes)) {
        if (node.dummy) {
            continue;
        }

        const newEdges = [];
        for (let neighbor of graph.edges[node.id]) {
            if (layers[neighbor] == layers[node.id] + 1) {
                newEdges.push(neighbor);
            } else {
                let prevVertex = node.id;
                for (let i = 1; i < layers[neighbor] - layers[node.id]; ++i) {
                    graph.edges[prevVertex].push(nextAvailableDummyId);
                    graph.nodes[nextAvailableDummyId] = {id: nextAvailableDummyId, dummy: true};
                    graph.edges[nextAvailableDummyId] = [];
                    layers[nextAvailableDummyId] = layers[node.id] + i;
                    prevVertex = nextAvailableDummyId;
                    nextAvailableDummyId -= 1;
                }
                graph.edges[prevVertex].push(neighbor);                
            }
        }
        graph.edges[node.id] = newEdges;
    }
}

function makeAdjacencyMatrix(prevLayer, nextLayer, graph) {
    const adjMatrix = {};
    for (let fromVertex of prevLayer) {
        adjMatrix[fromVertex] = {};
        for (let toVertex of nextLayer) {
            adjMatrix[fromVertex][toVertex] = 0;
        }
        for (let toVertex of graph.edges[fromVertex]) {
            adjMatrix[fromVertex][toVertex] = 1;
        }
    }
    return adjMatrix;
}

function transposeMatrix(matrix) {
    if (matrix.length == 0) {
        return matrix;
    }
    const rowNodes = Object.keys(matrix);
    const columnNodes = Object.keys(matrix[rowNodes[0]]);
    let transposedMatrix = {};
    for (let col of columnNodes) {
        transposedMatrix[col] = {};
        for (let row of rowNodes) {
            transposedMatrix[col][row] = matrix[row][col];
        }
    }
    return transposedMatrix;
}

function computeRowBarycenters(adjMatrix, rowNodes, columnNodes, weights) {
    if (typeof weights === 'undefined') {
        weights = [];
        for (let j = 0; j < columnNodes.length; ++j) {
            weights.push(j);
        }
    }

    const barycenters = [];
    for (let i = 0; i < rowNodes.length; ++i) {
        const fromVertex = rowNodes[i];
        let barycenter = 0;
        let total = 0;
        for (let j = 0; j < columnNodes.length; ++j) {
            const toVertex = columnNodes[j];
            barycenter += weights[j] * adjMatrix[fromVertex][toVertex];
            total += adjMatrix[fromVertex][toVertex];
        }
        if (total > 0) {
            barycenter /= total;
        }
        barycenters.push([barycenter, i]);
    }
    return barycenters;
}

function orderByBarycenter(adjMatrix, rowNodes, columnNodes) {
    const barycenters = computeRowBarycenters(adjMatrix, rowNodes, columnNodes);
    
    barycenters.sort(function(lhs, rhs) { 
        if (lhs[0] < rhs[0]) {
            return -1;
        }
        if (lhs[0] == rhs[0] && lhs[1] < rhs[1]) {
            return -1;
        }
        return 1;
    });

    const newOrder = []
    for (let entry of barycenters) {
        newOrder.push(rowNodes[entry[1]]);
    }
    return newOrder;
}

function tryUpdateBest(bestOrder, bestCrossings, nodesByLayer, graph) {
    let crossings = countCrossings(nodesByLayer, graph);
    if (crossings < bestCrossings) {
        bestOrder = nodesByLayer.slice();
        bestCrossings = crossings;
    }
    return [bestOrder, bestCrossings];
}

function countCrossingsVertexPair(adjMatrix, nextLayerOrder, firstVertex, secondVertex) {
    let crossings = 0;
    for (let i = 0; i < nextLayerOrder.length; ++i) {
        for (let j = i + 1; j < nextLayerOrder.length; ++j) {
            crossings += adjMatrix[firstVertex][nextLayerOrder[j]] * adjMatrix[secondVertex][nextLayerOrder[i]];
        }
    }
    return crossings;
}

function countCrossings(nodesByLayer, graph) {
    let crossings = 0;
    for (let layer = 1; layer < nodesByLayer.length; ++layer) {
        const adjMatrix = makeAdjacencyMatrix(nodesByLayer[layer - 1], nodesByLayer[layer], graph);

        for (let i = 0; i < nodesByLayer[layer - 1].length; ++i) {
            for (let j = i + 1; j < nodesByLayer[layer - 1].length; ++j) {
                crossings += countCrossingsVertexPair(adjMatrix, nodesByLayer[layer], nodesByLayer[layer - 1][i], nodesByLayer[layer - 1][j])
            }
        }
    }
    return crossings;
}

function phaseOne(bestOrder, bestCrossings, nodesByLayer, graph) {
    for (let iter = 0; iter < 5; ++iter) {
        for (let layer = 1; layer < nodesByLayer.length; ++layer) {
            let adjMatrix = makeAdjacencyMatrix(nodesByLayer[layer - 1], nodesByLayer[layer], graph);
            adjMatrix = transposeMatrix(adjMatrix);
            
            nodesByLayer[layer] = orderByBarycenter(adjMatrix, nodesByLayer[layer], nodesByLayer[layer - 1]);

            [bestOrder, bestCrossings] = tryUpdateBest(bestOrder, bestCrossings, nodesByLayer, graph);
        }

        for (let layer = nodesByLayer.length - 1; layer > 0; --layer) {
            let adjMatrix = makeAdjacencyMatrix(nodesByLayer[layer - 1], nodesByLayer[layer], graph);

            nodesByLayer[layer - 1] = orderByBarycenter(adjMatrix, nodesByLayer[layer - 1], nodesByLayer[layer]);
            
            [bestOrder, bestCrossings] = tryUpdateBest(bestOrder, bestCrossings, nodesByLayer, graph);
        }
    }
    return [bestOrder, bestCrossings, nodesByLayer];
}

function reverseSameValues(values) {
    if (values.length == 0) {
        return values;
    }

    let start = 0;
    let end = 1;
    while (end < values.length) {
        if (values[start][0] != values[end][0]) {
            values.splice(start, 0, ...values.splice(start, end - start).reverse());
            start = end;
        } else {
            end += 1;
        }
    }
    values.splice(start, 0, ...values.splice(start, end - start).reverse());
    return values;
}

function isStrictlyIncreasingSequence(values) {
    if (values.length == 0) {
        return true;
    }

    for (let i = 1; i < values.length; ++i) {
        if (values[i - 1][0] >= values[i][0]) {
            return false;
        }
    }
    return true;
}

function phaseTwo(bestOrder, bestCrossings, nodesByLayer, graph) {
    for (let iter = 0; iter < 5; ++iter) {
        for (let layer = 1; layer < nodesByLayer.length; ++layer) {
            let adjMatrix = makeAdjacencyMatrix(nodesByLayer[layer - 1], nodesByLayer[layer], graph);
            let rowBarycenters = computeRowBarycenters(adjMatrix, nodesByLayer[layer - 1], nodesByLayer[layer]);
            
            rowBarycenters = reverseSameValues(rowBarycenters);
            
            const newOrder = [];
            for (let entry of rowBarycenters) {
                newOrder.push(nodesByLayer[layer - 1][entry[1]]);
            }
            nodesByLayer[layer - 1] = newOrder;

            adjMatrix = makeAdjacencyMatrix(nodesByLayer[layer - 1], nodesByLayer[layer], graph);
            adjMatrix = transposeMatrix(adjMatrix);
            let columnBarycenters = computeRowBarycenters(adjMatrix, nodesByLayer[layer], nodesByLayer[layer - 1]);
            if (!isStrictlyIncreasingSequence(columnBarycenters)) {
                [bestOrder, bestCrossings, nodesByLayer] = phaseOne(bestOrder, bestCrossings, nodesByLayer, graph);
            }
        }

        for (let layer = nodesByLayer.length - 1; layer > 0; --layer) {
            let adjMatrix = makeAdjacencyMatrix(nodesByLayer[layer - 1], nodesByLayer[layer], graph);
            adjMatrix = transposeMatrix(adjMatrix);
            let columnBarycenters = computeRowBarycenters(adjMatrix, nodesByLayer[layer], nodesByLayer[layer - 1]);
            
            columnBarycenters = reverseSameValues(columnBarycenters);
            
            const newOrder = [];
            for (let entry of columnBarycenters) {
                newOrder.push(nodesByLayer[layer][entry[1]]);
            }
            nodesByLayer[layer] = newOrder;

            adjMatrix = makeAdjacencyMatrix(nodesByLayer[layer - 1], nodesByLayer[layer], graph);
            let rowBarycenters = computeRowBarycenters(adjMatrix, nodesByLayer[layer - 1], nodesByLayer[layer]);
            if (!isStrictlyIncreasingSequence(rowBarycenters)) {
                [bestOrder, bestCrossings, nodesByLayer] = phaseOne(bestOrder, bestCrossings, nodesByLayer, graph);
            }
        }
    }
    return [bestOrder, bestCrossings, nodesByLayer];
}

function edgeCrossingMinimization(graph, nodesByLayer) {
    let bestOrder = nodesByLayer.slice();
    let bestCrossings = countCrossings(bestOrder, graph);

    [bestOrder, bestCrossings, nodesByLayer] = phaseOne(bestOrder, bestCrossings, nodesByLayer, graph);
    
    [bestOrder, bestCrossings, nodesByLayer] = phaseTwo(bestOrder, bestCrossings, nodesByLayer, graph);
    return bestOrder;
}

function placementIteration(prevLayer, curLayer, graph, isDown) {
    const rowNodes = isDown ? prevLayer : curLayer;
    const columnNodes = isDown ? curLayer : prevLayer;

    let adjMatrix = makeAdjacencyMatrix(rowNodes, columnNodes, graph);
    if (!isDown) {
        adjMatrix = transposeMatrix(adjMatrix);
    }
    const priorities = [];
    for (let i = 0; i < curLayer.length; ++i) {
        const toVertex = curLayer[i];
        // if (graph.nodes[toVertex].dummy) {
        //     priorities.push([prevLayer.length + 1, i]);
        // } else {
            let degree = 0;
            for (let fromVertex of prevLayer) {
                degree += adjMatrix[fromVertex][toVertex];
            }
            priorities.push([degree, i]);
        //}
    }
    priorities.sort(function(lhs, rhs) { 
        if (lhs[0] > rhs[0]) {
            return -1;
        }
        if (lhs[0] == rhs[0] && lhs[1] < rhs[1]) {
            return -1;
        }
        return 1;
    });

    adjMatrix = transposeMatrix(adjMatrix);

    const weights = [];
    for (const vertex of prevLayer) {
        weights.push(graph.nodes[vertex].x);
    }
    const barycenters = computeRowBarycenters(adjMatrix, curLayer, prevLayer, weights);

    for (let entry of priorities) {
        const priority = entry[0];
        const curIndex = entry[1];
        const vertex = curLayer[curIndex];
        const barycenter = Math.round(barycenters[curIndex][0]);
        if (barycenter == graph.nodes[vertex].x) {
            continue;
        }

        if (barycenter < graph.nodes[vertex].x) {
            if (curIndex == 0) {
                graph.nodes[vertex].x = barycenter;
                continue;
            }

            let maxIndex = -1;
            for (let item of priorities) {
                if (item[0] >= priority && item[1] < curIndex) {
                    maxIndex = Math.max(maxIndex, item[1]);
                }
            }

            if (maxIndex == -1 || graph.nodes[curLayer[maxIndex]].x < barycenter - (curIndex - maxIndex) * HORIZONTAL_SCALE) {
                graph.nodes[vertex].x = barycenter;
                for (let l = 1; l <= curIndex; ++l) {
                    if (barycenter - l * HORIZONTAL_SCALE >= graph.nodes[curLayer[curIndex - l]].x) {
                        break;
                    }
                    graph.nodes[curLayer[curIndex - l]].x = barycenter - l * HORIZONTAL_SCALE;
                }
            } else {
                for (let l = 1; l <= curIndex - maxIndex; ++l) {
                    graph.nodes[curLayer[maxIndex + l]].x = graph.nodes[curLayer[maxIndex]].x + l * HORIZONTAL_SCALE;
                }
            }
        } else {
            if (curIndex == curLayer.length - 1) {
                graph.nodes[vertex].x = barycenter;
                continue;
            }

            let minIndex = curLayer.length;
            for (let item of priorities) {
                if (item[0] >= priority && item[1] > curIndex) {
                    minIndex = Math.min(minIndex, item[1]);
                }
            }
            if (minIndex == curLayer.length || graph.nodes[curLayer[minIndex]].x > barycenter + (minIndex - curIndex) * HORIZONTAL_SCALE) {
                graph.nodes[vertex].x = barycenter;
                for (let l = 1; l <= curLayer.length - curIndex - 1; ++l) {
                    if (barycenter + l * HORIZONTAL_SCALE <= graph.nodes[curLayer[curIndex + l]].x) {
                        break;
                    }
                    graph.nodes[curLayer[curIndex + l]].x = barycenter + l * HORIZONTAL_SCALE;
                }
            } else {
                for (let l = 1; l <= minIndex - curIndex; ++l) {
                    graph.nodes[curLayer[minIndex - l]].x = graph.nodes[curLayer[minIndex]].x - l * HORIZONTAL_SCALE;
                }
            }
        }
    }
} 

function computeHorizontalPlacement(graph, nodesByLayer) {
    for (let layer = 0; layer < nodesByLayer.length; ++layer) {
        for (let i = 0; i < nodesByLayer[layer].length; ++i) {
            graph.nodes[nodesByLayer[layer][i]].y = 50 + 100 * layer;
            graph.nodes[nodesByLayer[layer][i]].x = HORIZONTAL_SCALE * i;
        }
    }

    for (let layer = 1; layer < nodesByLayer.length; ++layer) {
        placementIteration(nodesByLayer[layer - 1], nodesByLayer[layer], graph, true);
    }

    for (let layer = nodesByLayer.length - 2; layer >= 0; --layer) {
        placementIteration(nodesByLayer[layer + 1], nodesByLayer[layer], graph, false);
    }

    for (let layer = 1; layer < nodesByLayer.length; ++layer) {
        placementIteration(nodesByLayer[layer - 1], nodesByLayer[layer], graph, true);
    }
}

function convertToEditorFormat(graph, withDummy) {
    let resultNodes = [];
    for (let node of Object.values(graph.nodes)) {
        if (withDummy || !node.dummy) {
            resultNodes.push(node);
        }
    }

    let resultLinks = [];
    for (let node of resultNodes) {
        for (let toVertex of graph.edges[node.id]) {
            let endVertex = toVertex;

            const interim = [];
            if (!withDummy) {
                while (graph.nodes[endVertex].dummy) {
                    interim.push(graph.nodes[endVertex]);
                    endVertex = graph.edges[endVertex][0];
                }
            }
            let source = node;
            let target = graph.nodes[endVertex];
            if (source.id < target.id) {
                resultLinks.push({source: source, target: target, left: false, right: true, interim: interim});
            } else {
                resultLinks.push({source: target, target: source, left: true, right: false, interim: interim});
            }
        }
    }

    return [resultNodes, resultLinks];
}