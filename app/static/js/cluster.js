fetch('js/data.json', {mode: 'no-cors'})
  .then(function(res) {
    return res.json()
  })
  .then(function(data) {
  // Options for the cytoscape.js visualization; we are using the 'cise' layout 
  let options = {
    name: 'cise',
    clusters: data.clusterInfo,
    animate: false,
    refresh: 10, 
    animationDuration: undefined,
    animationEasing: undefined,
    fit: true,
    padding: 20,
    nodeSeparation: 20,
    idealInterClusterEdgeLengthCoefficient: 1.4,
    allowNodesInsideCircle: false,
    maxRatioOfNodesInsideCircle: 0.1,
    springCoeff: 0.45,
    nodeRepulsion: 1,
    gravity: 0.25,
    gravityRange: 3.8, 
    ready: function(){}, 
    stop: function(){},
  };

  let rgb1 = data['color1'] 
  let rgb2 = data['color2']

  let val_range = data['highest_val'] - data['lowest_val']

  function getRgbForVal(d, a1, a2) {
    // finds an rgb d in between a1 and a2,
    return a1.map((p, i) => Math.floor(a1[i] + d * (a2[i] - a1[i])))
  } 

  let rgbToHex = function (rgb) {
		let hex = Number(rgb).toString(16);
		if (hex.length < 2) {
			hex = "0" + hex;
		}
		return hex;
	};

  let fullColorHex = function(r,g,b) {
		var red = rgbToHex(r);
		var green = rgbToHex(g);
		var blue = rgbToHex(b);
		return red+green+blue;
  };
  
  var cy = window.cy = cytoscape({
    container: document.getElementById('cy'),
    boxSelectionEnabled: false,
    autounselectify: true,

    style: [
      {
        selector: 'node',
        style: {
          'height': 6,
          'width': 6,
          'label': 'data(label)',
          'text-wrap': 'wrap',
          'font-size': '3px',
          'background-color' : function(ele)  {
            //for our color gradient, we want to find the range of property values, 
            //then see where the current property value falls in that range and color it accordingly
            let rgbColors = getRgbForVal((ele.data('prop_val') - data['lowest_val'])/val_range, rgb1, rgb2);
            let hexColor = '#' + fullColorHex(rgbColors[0], rgbColors[1], rgbColors[2]);
            return hexColor;
          },
          'shape' : function(ele) {
            // centroids are represented as stars
            if(ele.data('centroid')){
              return 'star';
            }
            // if the chemical has been reclustered, it is a triangle
            if(ele.data('reclustered')){
              return 'triangle'
            }
            return 'ellipse';
          }
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

  function makePopper(ele) {
    let ref = ele.popperRef(); // used only for positioning
  }

  cy.ready(function() {
    cy.elements().forEach(function(ele) {
      if(ele.data('type')== 'node') {
        makePopper(ele);
      }
    });
  });

  // Draw the current smile structure that the mouse is hovering over
  cy.elements().unbind("mouseover");
  cy.elements().bind("mouseover", event => {
    const smile = event.target.id()
    let smilesDrawer = new SmilesDrawer.Drawer({width: 400, height:400});
    SmilesDrawer.parse(smile, function(tree) {
      smilesDrawer.draw(tree, "drawing", "light", false);
    });
  });

  cy.elements().unbind("mouseout");

  // actually run and create the cytoscape visualization
  var layout = cy.layout( options );
  layout.run();

  // Download image button
  const img_options =
  {
    "bg": "white",
    "full": false,
    "maxWidth": 15000,
    "maxHeight": 11250,
  };

  var png64 = cy.png(img_options);
  
  $('<div class=\'text-center\'><a id="png" download>Download Image!</a></div>').insertBefore('#cy');
    $('#png').attr('href', png64);

  // create color squares for legend
  const color1 = '#' + fullColorHex(rgb1[0], rgb1[1], rgb1[2]);
  const color2 = '#' + fullColorHex(rgb2[0], rgb2[1], rgb2[2]);

  const sq1 = '<svg width="20" height="20"> <rect width="20" height="20" style="fill:' + color1 +';stroke-width:3;stroke:rgb(0,0,0)" /> </svg>';
  const sq2 = '<svg width="20" height="20"> <rect width="20" height="20" style="fill:' + color2 +';stroke-width:3;stroke:rgb(0,0,0)" /> </svg>';

  document.getElementById('square1').innerHTML = sq1;
  document.getElementById('square2').innerHTML = sq2;

  var resetButton = document.createElement("button");
  resetButton.innerHTML = "<div style=\'text-align: center\'>Reset Cluster View</div>";

  // Add a reset cluster view button
  resetButton.addEventListener ("click", function() {
    layout.stop()
    layout.run();
  });

  var legend = document.getElementById("legend");
  var legendParent = legend.parentNode;

  legendParent.insertBefore(resetButton, legend);

});
