const bamtaxSelectDataset = window.selectDataset;

function hideKronaDatasetControls() {
  const dropdown = document.getElementById('datasets');
  const controls = dropdown && dropdown.closest('div');
  if (controls) {
    controls.hidden = true;
  }
}

if (typeof bamtaxSelectDataset === 'function') {
  window.selectDataset = index => {
    hideKronaDatasetControls();
    bamtaxSelectDataset(index);
    window.parent.postMessage({
      type: 'bamtax-dataset-selected',
      index,
    }, '*');
  };
}

window.addEventListener('load', hideKronaDatasetControls);

window.addEventListener('message', event => {
  const message = event.data;
  if (
    event.source === window.parent &&
    message &&
    message.type === 'bamtax-dataset' &&
    Number.isInteger(message.index) &&
    typeof window.selectDataset === 'function'
  ) {
    window.selectDataset(message.index);
  }
});
