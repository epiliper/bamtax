const selector = document.getElementById('dataset');
const panels = [...document.querySelectorAll('.panel')];
const krona = document.querySelector('iframe');

function activePanel() {
  return document.getElementById(selector.value);
}

function syncKrona() {
  if (krona) {
    krona.contentWindow.postMessage({
      type: 'bamtax-dataset',
      index: Number(selector.value.slice(6)),
    }, '*');
  }
}

function selectDataset(index, notifyKrona) {
  const value = `panel-${index}`;
  if (!document.getElementById(value)) {
    return;
  }

  selector.value = value;
  panels.forEach(panel => panel.hidden = panel.id !== value);
  resetSearches();
  renderSankey(index);
  if (notifyKrona) {
    syncKrona();
  }
}

const sankeyView = document.getElementById('sankey-view');
const kronaView = document.getElementById('krona-view');
const showSankey = document.getElementById('show-sankey');
const showKrona = document.getElementById('show-krona');
const sankeySampleName = document.getElementById('sankey-sample-name');

function selectVisualization(sankey) {
  sankeyView.hidden = !sankey;
  kronaView.hidden = sankey;
  showSankey.classList.toggle('active', sankey);
  showKrona.classList.toggle('active', !sankey);
  showSankey.setAttribute('aria-pressed', String(sankey));
  showKrona.setAttribute('aria-pressed', String(!sankey));
  if (sankey) {
    scheduleSankeyRender();
  }
}

showSankey.addEventListener('click', () => selectVisualization(true));
showKrona.addEventListener('click', () => selectVisualization(false));

function renderSankey(index) {
  const svg = d3.select('#sankey');
  svg.selectAll('*').remove();
  hideSankeyTooltip();
  const dataset = sankeyDatasets[index];
  sankeySampleName.textContent = dataset?.name || '';
  if (!dataset || dataset.species.length === 0) {
    svg.attr('viewBox', '0 0 800 120').attr('width', 800).attr('height', 120)
      .append('text').attr('x', 400).attr('y', 65)
      .attr('text-anchor', 'middle').attr('fill', 'var(--muted)').text('No species calls passed the report filters.');
    return;
  }

  const { graph, width, height, labelX } = layoutSankey(dataset, Math.max(0, sankeyView.clientWidth - 2));
  svg.attr('viewBox', `0 0 ${width} ${height}`).attr('width', width).attr('height', height);
  svg.append('g').attr('fill', 'none').attr('stroke-opacity', .26).selectAll('path').data(graph.links).join('path')
    .attr('d', d3.sankeyLinkHorizontal()).attr('stroke', 'var(--sankey-link)').attr('stroke-width', link => Math.max(1, link.width));
  const speciesNodes = graph.nodes.filter(item => item.species).sort((a, b) => nodeCenter(a) - nodeCenter(b));
  const labelPositions = spreadLabels(speciesNodes.map(nodeCenter), 20, height - 20, 18);
  svg.append('g').selectAll('path').data(speciesNodes).join('path').attr('class', 'species-leader')
    .attr('d', (item, itemIndex) => {
      const startY = nodeCenter(item);
      const endY = labelPositions[itemIndex];
      const elbowX = item.x1 + 12;
      return `M${item.x1},${startY}H${elbowX}V${endY}H${labelX - 7}`;
    });
  const node = svg.append('g').selectAll('g').data(graph.nodes).join('g')
    .attr('class', item => `node${!item.species ? ' interactive-node' : ''}`)
    .attr('tabindex', item => item.species ? null : 0)
    .attr('role', 'img')
    .attr('aria-label', nodeDescription);
  node.append('rect').attr('x', item => item.x0).attr('y', item => item.y0)
    .attr('width', item => item.x1 - item.x0).attr('height', item => Math.max(1, item.y1 - item.y0))
    .attr('fill', item => item.species ? 'var(--sankey-species)' : 'var(--sankey-node)');
  const groupLabels = node.filter(item => !item.species).append('text')
    .attr('x', item => item.x0 - 10).attr('y', item => nodeCenter(item) - 3).attr('text-anchor', 'end');
  groupLabels.append('tspan').text(item => item.name);
  groupLabels.append('tspan').attr('x', item => item.x0 - 10).attr('dy', '1.2em').text(item => item.value);
  node.filter(item => !item.species)
    .on('focus', showSankeyTooltip)
    .on('blur', hideTransientSankeyTooltip)
    .on('click', toggleSankeyTooltip);
  svg.append('g').selectAll('text').data(speciesNodes).join('text')
    .attr('x', labelX).attr('y', (_, itemIndex) => labelPositions[itemIndex]).attr('dy', '.35em')
    .text(item => `${item.name} (${item.value})`);
  svg.on('click', hideSankeyTooltip);
}

function layoutSankey(dataset, availableWidth) {
  const groupOrder = [2, 2157, 2759, 10239, 0];
  const branchIds = [...new Set(dataset.species.map(species => species.lineage[0].taxId))]
    .sort((a, b) => groupOrder.indexOf(a) - groupOrder.indexOf(b));
  const branchGap = 28;
  const branchExtra = Math.max(0, branchIds.length - 1) * branchGap;
  const longestLabel = Math.max(...dataset.species.map(species => {
    const taxon = species.lineage[species.lineage.length - 1];
    return `${taxon.name} (${species.value})`.length;
  }));
  const labelWidth = Math.min(420, Math.max(210, longestLabel * 7 + 20));
  const width = Math.max(680, availableWidth, 2 * 125 + labelWidth + 140);
  const height = Math.max(300, dataset.species.length * 34 + 40 + branchExtra);
  const labelX = width - labelWidth;
  const plotRight = labelX - 32;
  const nodesById = new Map();
  const linksById = new Map();
  dataset.species.forEach(species => {
    let source = null;
    const branchIndex = branchIds.indexOf(species.lineage[0].taxId);
    species.lineage.forEach((taxon, lineageIndex) => {
      const target = `${lineageIndex === 0 ? 'group' : 'taxon'}-${taxon.taxId}`;
      if (!nodesById.has(target)) {
        nodesById.set(target, {
          id: target,
          name: taxon.name,
          taxId: taxon.taxId,
          rank: taxon.rank,
          branchIndex,
          species: lineageIndex === species.lineage.length - 1,
          fixedValue: 0,
        });
      } else if (lineageIndex === species.lineage.length - 1) {
        nodesById.get(target).species = true;
      }
      nodesById.get(target).fixedValue += species.value;
      if (source !== null) {
        const linkId = `${source}>${target}`;
        const link = linksById.get(linkId);
        if (link) {
          link.value += species.value;
        } else {
          linksById.set(linkId, { source, target, value: species.value });
        }
      }
      source = target;
    });
  });
  const nodes = [...nodesById.values()];
  const links = [...linksById.values()];
  const sankey = d3.sankey().nodeId(node => node.id).nodeWidth(24).nodePadding(18)
    .nodeSort((a, b) => a.branchIndex - b.branchIndex || b.fixedValue - a.fixedValue || a.name.localeCompare(b.name) || a.taxId - b.taxId)
    .extent([[140, 20], [plotRight, height - 20 - branchExtra]]);
  const graph = sankey({ nodes, links });
  graph.nodes.forEach(item => {
    const offset = item.branchIndex * branchGap;
    item.y0 += offset;
    item.y1 += offset;
  });
  // D3's relaxation can leave attachment order out of sync with final node positions.
  graph.nodes.forEach(item => {
    item.sourceLinks.sort((a, b) => a.target.y0 - b.target.y0 || a.index - b.index);
    item.targetLinks.sort((a, b) => a.source.y0 - b.source.y0 || a.index - b.index);
  });
  sankey.update(graph);
  return { graph, width, height, labelX };
}

function nodeCenter(node) {
  return (node.y0 + node.y1) / 2;
}

function nodeDescription(node) {
  return `${node.name} (${node.rank || 'unranked'}, taxid ${node.taxId}): ${node.value}`;
}

function spreadLabels(positions, top, bottom, gap) {
  const labels = positions.map(position => Math.max(top, position));
  for (let index = 1; index < labels.length; index += 1) {
    labels[index] = Math.max(labels[index], labels[index - 1] + gap);
  }
  if (labels.length && labels[labels.length - 1] > bottom) {
    labels[labels.length - 1] = bottom;
    for (let index = labels.length - 2; index >= 0; index -= 1) {
      labels[index] = Math.min(labels[index], labels[index + 1] - gap);
    }
  }
  return labels;
}

const sankeyTooltip = document.getElementById('sankey-tooltip');
let pinnedSankeyNode = null;

function showSankeyTooltip(event, node) {
  if (pinnedSankeyNode && pinnedSankeyNode !== event.currentTarget) {
    return;
  }
  sankeyTooltip.textContent = nodeDescription(node);
  sankeyTooltip.hidden = false;
  moveSankeyTooltip(event);
}

function moveSankeyTooltip(event) {
  if (pinnedSankeyNode && pinnedSankeyNode !== event.currentTarget) {
    return;
  }
  const viewBounds = sankeyView.getBoundingClientRect();
  const targetBounds = event.currentTarget.getBoundingClientRect();
  const clientX = event.clientX || targetBounds.right;
  const clientY = event.clientY || targetBounds.top + targetBounds.height / 2;
  const left = clientX - viewBounds.left + sankeyView.scrollLeft + 12;
  const top = clientY - viewBounds.top + sankeyView.scrollTop + 12;
  sankeyTooltip.style.left = `${Math.min(left, sankeyView.scrollLeft + sankeyView.clientWidth - sankeyTooltip.offsetWidth - 8)}px`;
  sankeyTooltip.style.top = `${top}px`;
}

function hideSankeyTooltip() {
  if (pinnedSankeyNode) {
    pinnedSankeyNode.classList.remove('pinned');
    pinnedSankeyNode = null;
  }
  sankeyTooltip.hidden = true;
}

function hideTransientSankeyTooltip() {
  if (!pinnedSankeyNode) {
    sankeyTooltip.hidden = true;
  }
}

function toggleSankeyTooltip(event, node) {
  event.stopPropagation();
  if (pinnedSankeyNode === event.currentTarget) {
    hideSankeyTooltip();
    return;
  }
  if (pinnedSankeyNode) {
    pinnedSankeyNode.classList.remove('pinned');
  }
  pinnedSankeyNode = event.currentTarget;
  pinnedSankeyNode.classList.add('pinned');
  sankeyTooltip.textContent = nodeDescription(node);
  sankeyTooltip.hidden = false;
  moveSankeyTooltip(event);
}

document.addEventListener('keydown', event => {
  if (event.key === 'Escape') {
    hideSankeyTooltip();
    hideSearchResults();
  }
});

let sankeyRenderFrame;
function scheduleSankeyRender() {
  cancelAnimationFrame(sankeyRenderFrame);
  sankeyRenderFrame = requestAnimationFrame(() => {
    if (!sankeyView.hidden) {
      renderSankey(Number(selector.value.slice(6)));
    }
  });
}

if ('ResizeObserver' in window) {
  let previousWidth = Math.round(sankeyView.getBoundingClientRect().width);
  new ResizeObserver(entries => {
    const width = Math.round(entries[0].contentRect.width);
    if (width !== previousWidth) {
      previousWidth = width;
      scheduleSankeyRender();
    }
  }).observe(sankeyView);
} else {
  window.addEventListener('resize', scheduleSankeyRender);
}

let sankeyPan = null;
sankeyView.addEventListener('pointerdown', event => {
  if (event.pointerType !== 'mouse' || event.button !== 0 || event.target.closest('.interactive-node')) {
    return;
  }
  sankeyPan = {
    pointerId: event.pointerId,
    x: event.clientX,
    y: event.clientY,
    left: sankeyView.scrollLeft,
    top: sankeyView.scrollTop,
  };
  sankeyView.setPointerCapture(event.pointerId);
  sankeyView.classList.add('panning');
});

sankeyView.addEventListener('pointermove', event => {
  if (!sankeyPan || event.pointerId !== sankeyPan.pointerId) {
    return;
  }
  sankeyView.scrollLeft = sankeyPan.left - (event.clientX - sankeyPan.x);
  sankeyView.scrollTop = sankeyPan.top - (event.clientY - sankeyPan.y);
});

function stopSankeyPan(event) {
  if (sankeyPan && event.pointerId === sankeyPan.pointerId) {
    sankeyPan = null;
    sankeyView.classList.remove('panning');
  }
}

sankeyView.addEventListener('pointerup', stopSankeyPan);
sankeyView.addEventListener('pointercancel', stopSankeyPan);

selector.addEventListener('change', () => selectDataset(Number(selector.value.slice(6)), true));

window.addEventListener('message', event => {
  const message = event.data;
  if (!krona || event.source !== krona.contentWindow || !message) {
    return;
  }
  if (message.type === 'bamtax-dataset-selected' && Number.isInteger(message.index)) {
    selectDataset(message.index, false);
  }
});

if (krona) {
  krona.addEventListener('load', syncKrona);
}

const reportTabs = [...document.querySelectorAll('.report-tabs [role="tab"]')];
let heatmapRendered = false;
let heatmapRenderFrame;
function scheduleHeatmapRender() {
  cancelAnimationFrame(heatmapRenderFrame);
  heatmapRenderFrame = requestAnimationFrame(() => {
    if (heatmapRendered && !document.getElementById('heatmap-panel').hidden) {
      renderHeatmap(sankeyDatasets);
    }
  });
}
if ('ResizeObserver' in window) {
  let previousWidth = 0;
  new ResizeObserver(entries => {
    const width = Math.round(entries[0].contentRect.width);
    if (width !== previousWidth) {
      previousWidth = width;
      scheduleHeatmapRender();
    }
  }).observe(document.getElementById('heatmap-view'));
} else {
  window.addEventListener('resize', scheduleHeatmapRender);
}
function selectReportTab(tab) {
  reportTabs.forEach(item => {
    const selected = item === tab;
    item.classList.toggle('active', selected);
    item.setAttribute('aria-selected', String(selected));
    item.tabIndex = selected ? 0 : -1;
    document.getElementById(item.getAttribute('aria-controls')).hidden = !selected;
  });
  hideSankeyTooltip();
  document.getElementById('heatmap-tooltip').hidden = true;
  if (tab.id === 'tab-heatmap') {
    renderHeatmap(sankeyDatasets);
    heatmapRendered = true;
  } else if (tab.id === 'tab-explorer' && !sankeyView.hidden) {
    scheduleSankeyRender();
  }
}
reportTabs.forEach((tab, index) => {
  tab.addEventListener('click', () => selectReportTab(tab));
  tab.addEventListener('keydown', event => {
    if (!['ArrowLeft', 'ArrowRight', 'Home', 'End'].includes(event.key)) return;
    event.preventDefault();
    const next = event.key === 'Home' ? 0 : event.key === 'End' ? reportTabs.length - 1
      : (index + (event.key === 'ArrowRight' ? 1 : -1) + reportTabs.length) % reportTabs.length;
    selectReportTab(reportTabs[next]);
    reportTabs[next].focus();
  });
});

renderSankey(0);

document.getElementById('expand').addEventListener('click', () => {
  activePanel().querySelectorAll('details').forEach(node => node.open = true);
});

document.getElementById('collapse').addEventListener('click', () => {
  activePanel().querySelectorAll('details').forEach(node => node.open = false);
});

const extractDialog = document.getElementById('extract-dialog');
const extractCommand = document.getElementById('extract-command');
const copyStatus = document.getElementById('copy-status');

async function copyExtractCommand() {
  let copied = false;
  if (navigator.clipboard) {
    try {
      await navigator.clipboard.writeText(extractCommand.value);
      copied = true;
    } catch (_) {}
  }

  if (!copied) {
    extractCommand.select();
    copied = document.execCommand('copy');
  }

  copyStatus.textContent = copied ? 'Copied to clipboard' : 'Select and copy the command manually';
}

document.querySelectorAll('.extract-button').forEach(button => {
  button.addEventListener('click', event => {
    event.preventDefault();
    event.stopPropagation();
    extractCommand.value = button.dataset.command;
    copyStatus.textContent = '';

    if (extractDialog.showModal) {
      extractDialog.showModal();
    } else {
      extractDialog.setAttribute('open', '');
    }

    copyExtractCommand();
  });
});

document.getElementById('copy-command').addEventListener('click', copyExtractCommand);

const referencesDialog = document.getElementById('references-dialog');
const referencesTitle = document.getElementById('references-title');
const referencesList = document.getElementById('references-list');

function openReferencesDialog(taxon, references) {
  referencesTitle.textContent = `References for ${taxon}`;
  referencesList.replaceChildren();

  references.split('\n').filter(Boolean).forEach(reference => {
    const item = document.createElement('li');
    item.textContent = reference;
    referencesList.appendChild(item);
  });

  if (referencesDialog.showModal) {
    referencesDialog.showModal();
  } else {
    referencesDialog.setAttribute('open', '');
  }
}

document.querySelectorAll('.references-button').forEach(button => {
  button.addEventListener('click', event => {
    event.preventDefault();
    event.stopPropagation();
    openReferencesDialog(button.dataset.taxon, button.dataset.references);
  });
});

function activeTaxonEntries() {
  return [...activePanel().querySelectorAll('.search-taxon')];
}

const taxonSearch = document.getElementById('taxon-search');
const referenceSearch = document.getElementById('reference-search');
const taxonSearchResults = document.getElementById('taxon-search-results');
const referenceSearchResults = document.getElementById('reference-search-results');

function addSearchResult(results, label, activate) {
  const result = document.createElement('button');
  result.type = 'button';
  result.textContent = label;
  result.addEventListener('click', () => {
    hideSearchResults();
    activate();
  });
  results.appendChild(result);
}

function focusTaxonRow(entry) {
  let ancestor = entry.parentElement;
  while (ancestor) {
    if (ancestor.tagName === 'DETAILS') {
      ancestor.open = true;
    }
    ancestor = ancestor.parentElement;
  }

  const row = entry.closest('.row');
  if (!row) {
    return;
  }
  document.querySelectorAll('.row.search-focus').forEach(item => item.classList.remove('search-focus'));
  row.classList.add('search-focus');
  row.setAttribute('tabindex', '-1');
  row.focus({ preventScroll: true });
  row.scrollIntoView({ behavior: 'smooth', block: 'center', inline: 'nearest' });
}

function renderTaxonSearch() {
  const query = taxonSearch.value.trim().toLocaleLowerCase();
  taxonSearchResults.replaceChildren();
  if (!query) {
    taxonSearchResults.hidden = true;
    return;
  }
  const matches = activeTaxonEntries().filter(entry =>
    entry.dataset.taxon.toLocaleLowerCase().includes(query) || entry.dataset.taxId.includes(query)
  ).slice(0, 50);
  matches.forEach(entry => addSearchResult(
    taxonSearchResults,
    `${entry.dataset.taxon} — ${entry.dataset.rank || 'unranked'}`,
    () => focusTaxonRow(entry),
  ));
  showSearchResultState(taxonSearchResults, matches.length);
}

function renderReferenceSearch() {
  const query = referenceSearch.value.trim().toLocaleLowerCase();
  referenceSearchResults.replaceChildren();
  if (!query) {
    referenceSearchResults.hidden = true;
    return;
  }
  const matches = [];
  activeTaxonEntries().forEach(entry => {
    entry.dataset.references.split('\n').filter(reference => reference.toLocaleLowerCase().includes(query))
      .forEach(reference => matches.push({ reference, entry }));
  });
  matches.slice(0, 50).forEach(match => addSearchResult(
    referenceSearchResults,
    `${match.reference} — ${match.entry.dataset.taxon}`,
    () => openReferencesDialog(match.entry.dataset.taxon, match.entry.dataset.references),
  ));
  showSearchResultState(referenceSearchResults, matches.length);
}

function showSearchResultState(results, count) {
  if (count === 0) {
    const empty = document.createElement('span');
    empty.className = 'search-empty';
    empty.textContent = 'No matches in this sample';
    results.appendChild(empty);
  }
  results.hidden = false;
}

function hideSearchResults() {
  taxonSearchResults.hidden = true;
  referenceSearchResults.hidden = true;
}

function resetSearches() {
  taxonSearch.value = '';
  referenceSearch.value = '';
  taxonSearchResults.replaceChildren();
  referenceSearchResults.replaceChildren();
  hideSearchResults();
}

taxonSearch.addEventListener('input', renderTaxonSearch);
taxonSearch.addEventListener('focus', renderTaxonSearch);
referenceSearch.addEventListener('input', renderReferenceSearch);
referenceSearch.addEventListener('focus', renderReferenceSearch);
document.addEventListener('click', event => {
  if (!event.target.closest('.search-field')) {
    hideSearchResults();
  }
});
