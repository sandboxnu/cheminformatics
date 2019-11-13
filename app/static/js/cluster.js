fetch('js/data.json', {mode: 'no-cors'})
  .then(function(res) {
    return res.json()
  })
  .then(function(data) {
  let options = {
    name: 'circle',

    fit: true, 
    padding: 30,
    boundingBox: undefined, 
    avoidOverlap: true, 
    nodeDimensionsIncludeLabels: false,
    radius: undefined, 
    startAngle: 3 / 2 * Math.PI, 
    sweep: undefined,
    clockwise: true, 
    sort: undefined, 
    animationEasing: undefined, 
    animateFilter: function ( node, i ){ return true; }, 
    ready: undefined,
    stop: undefined, 
    transform: function (node, position ){ return position; } 

  };

  var cy = window.cy = cytoscape({
    container: document.getElementById('cy'),

    boxSelectionEnabled: false,
    autounselectify: true,

    style: [
      {
        selector: 'node',
        style: {
          'height': 20,
          'width': 20,
          'label': 'data(label)',
        }
      },

      {
        selector: 'edge',
        style: {
          'curve-style': 'haystack',
          'haystack-radius': 0,
          'width': 5,
          'opacity': 0.5,
          'line-color': '#f2f08c'
        }
      }
    ],

    elements: data
  });

  var layout = cy.layout( options );

  layout.run();
});