fetch('js/data.json', {mode: 'no-cors'})
  .then(function(res) {
    return res.json()
  })
  .then(function(data) {
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



  let rgb1 = data['color1']   // red
  let rgb2 = data['color2']  // yellow

  let val_range = data['highest_val'] - data['lowest_val']
  function getPoint(d, a1, a2) {
  // find a color d% between a1 and a2
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
            let colors = getPoint((ele.data('prop_val') - data['lowest_val'])/val_range, rgb1, rgb2);
            let rgbColor = '#' + fullColorHex(colors[0], colors[1], colors[2]);
            return rgbColor;
          },
          'shape' : function(ele) {
            if(ele.data('centroid')){
              return 'star';
            }
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

  cy.elements().unbind("mouseover");
  cy.elements().bind("mouseover", event => {
    const smile = event.target.id()
    let smilesDrawer = new SmilesDrawer.Drawer({width: 400, height:400});
    SmilesDrawer.parse(smile, function(tree) {
      // Draw to the canvas
      smilesDrawer.draw(tree, "drawing", "light", false);
      // Alternatively, draw to SVG:
      // svgDrawer.draw(tree, 'output-svg', 'dark', false);
    });
  });

  cy.elements().unbind("mouseout");
  cy.elements().bind("mouseout", event => {
    //event.target.tippy.hide()
  });

  var layout = cy.layout( options );
  layout.run();
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

  
  var resetButton = document.createElement("button");
  resetButton.innerHTML = "<div style=\'text-align: center\'>Reset Cluster View</div>";
 
  resetButton.addEventListener ("click", function() {
    layout.stop()
    layout.run();
  });

  var legend = document.getElementById("legend");
  var legendParent = legend.parentNode;

  legendParent.insertBefore(resetButton, legend);

});
