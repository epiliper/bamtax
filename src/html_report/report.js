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
  renderSankey(index);
  if (notifyKrona) {
    syncKrona();
  }
}

const sankeyView = document.getElementById('sankey-view');
const kronaView = document.getElementById('krona-view');
const showSankey = document.getElementById('show-sankey');
const showKrona = document.getElementById('show-krona');

function selectVisualization(sankey) {
  sankeyView.hidden = !sankey;
  kronaView.hidden = sankey;
  showSankey.classList.toggle('active', sankey);
  showKrona.classList.toggle('active', !sankey);
  showSankey.setAttribute('aria-pressed', String(sankey));
  showKrona.setAttribute('aria-pressed', String(!sankey));
}

showSankey.addEventListener('click', () => selectVisualization(true));
showKrona.addEventListener('click', () => selectVisualization(false));

function renderSankey(index) {
  const svg = d3.select('#sankey');
  svg.selectAll('*').remove();
  const dataset = sankeyDatasets[index];
  if (!dataset || dataset.species.length === 0) {
    svg.attr('viewBox', '0 0 800 120').append('text').attr('x', 400).attr('y', 65)
      .attr('text-anchor', 'middle').attr('fill', '#65716c').text('No species calls passed the report filters.');
    return;
  }

  const width = 900;
  const height = Math.max(300, dataset.species.length * 32);
  const sourceId = 'sample';
  const nodes = [{ id: sourceId, name: dataset.name, source: true }]
    .concat(dataset.species.map(species => ({ id: `taxon-${species.taxId}`, name: species.name, taxId: species.taxId })));
  const links = dataset.species.map(species => ({ source: sourceId, target: `taxon-${species.taxId}`, value: species.value }));
  const graph = d3.sankey().nodeId(node => node.id).nodeWidth(18).nodePadding(12)
    .extent([[150, 20], [width - 200, height - 20]])({ nodes, links });
  svg.attr('viewBox', `0 0 ${width} ${height}`);
  svg.append('g').attr('fill', 'none').attr('stroke-opacity', .4).selectAll('path').data(graph.links).join('path')
    .attr('d', d3.sankeyLinkHorizontal()).attr('stroke', '#176b52').attr('stroke-width', link => Math.max(1, link.width))
    .append('title').text(link => `${link.target.name}: ${link.value}`);
  const node = svg.append('g').selectAll('g').data(graph.nodes).join('g');
  node.append('rect').attr('x', item => item.x0).attr('y', item => item.y0)
    .attr('width', item => item.x1 - item.x0).attr('height', item => Math.max(1, item.y1 - item.y0))
    .attr('fill', item => item.source ? '#176b52' : '#efb55f');
  node.append('text').attr('x', item => item.source ? item.x0 - 8 : item.x1 + 8)
    .attr('y', item => (item.y0 + item.y1) / 2).attr('dy', '.35em')
    .attr('text-anchor', item => item.source ? 'end' : 'start')
    .text(item => item.source ? item.name : `${item.name} (${item.value})`);
  node.append('title').text(item => item.source ? item.name : `${item.name} (taxid ${item.taxId}): ${item.value}`);
}

selector.addEventListener('change', () => selectDataset(Number(selector.value.slice(6)), true));

window.addEventListener('message', event => {
  const message = event.data;
  if (
    krona &&
    event.source === krona.contentWindow &&
    message &&
    message.type === 'bamtax-dataset-selected' &&
    Number.isInteger(message.index)
  ) {
    selectDataset(message.index, false);
  }
});

if (krona) {
  krona.addEventListener('load', syncKrona);
}

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

document.querySelectorAll('.references-button').forEach(button => {
  button.addEventListener('click', event => {
    event.preventDefault();
    event.stopPropagation();
    referencesTitle.textContent = `References for ${button.dataset.taxon}`;
    referencesList.replaceChildren();

    button.dataset.references.split('\n').filter(Boolean).forEach(reference => {
      const item = document.createElement('li');
      item.textContent = reference;
      referencesList.appendChild(item);
    });

    if (referencesDialog.showModal) {
      referencesDialog.showModal();
    } else {
      referencesDialog.setAttribute('open', '');
    }
  });
});
