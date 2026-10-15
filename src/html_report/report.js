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

selector.addEventListener('change', () => {
  panels.forEach(panel => panel.hidden = panel.id !== selector.value);
  syncKrona();
});

if (krona) {
  krona.addEventListener('load', syncKrona);
}

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
