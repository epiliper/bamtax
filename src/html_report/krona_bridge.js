window.addEventListener('message', event => {
  const message = event.data;
  if (
    event.source === window.parent &&
    message &&
    message.type === 'bamtax-dataset' &&
    Number.isInteger(message.index) &&
    typeof selectDataset === 'function'
  ) {
    selectDataset(message.index);
  }
});
