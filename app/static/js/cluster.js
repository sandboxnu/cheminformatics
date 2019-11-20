fetch('js/data.json', {mode: 'no-cors'})
  .then(function(res) {
    return res.json()
  })
  .then(function(data) {
  let options = {
    name: 'cise',

    clusters: function(node) { return 1},
    animate: false,
    
    refresh: 10, 
    animationDuration: undefined,
    animationEasing: undefined,
    fit: true,
    padding: 30,
    nodeSeparation: 20,
    idealInterClusterEdgeLengthCoefficient: 1.4,
    allowNodesInsideCircle: false,
    maxRatioOfNodesInsideCircle: 0.1,
    springCoeff: 0.45,
    nodeRepulsion: 4500,
    gravity: 0.25,
    gravityRange: 3.8, 
    ready: function(){}, 
    stop: function(){},
  };

  var cy = window.cy = cytoscape({
    container: document.getElementById('cy'),

    boxSelectionEnabled: false,
    autounselectify: true,

    style: [
      {
        selector: 'node',
        style: {
          'height': 10,
          'width': 10,
          'label': 'data(label)',
          'text-wrap': 'wrap'
        }
      },

      {
        selector: 'edge',
        style: {
          'curve-style': 'haystack',
          'haystack-radius': 0,
          'width': 2,
          'opacity': 0.5,
          'line-color': '#f2f08c'
        }
      }
    ],

    elements: data
  });

  var layout = cy.layout( options );

  layout.run();

  var png64 = cy.png();

  $('<div class=\'text-center\'><a id="png" download>Download Image!</a></div>').insertBefore('#cy');
    $('#png').attr('href', png64);
});
