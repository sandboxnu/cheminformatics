fetch('js/data.json', {mode: 'no-cors'})
  .then(function(res) {
    return res.json()
  })
  .then(function(data){
    console.log(data);
    let rgb1 = data['color1']   // red
    let rgb2 = data['color2']  // yellow'

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
    
    const color1 = '#' + fullColorHex(rgb1[0], rgb1[1], rgb1[2]);
    const color2 = '#' + fullColorHex(rgb2[0], rgb2[1], rgb2[2]);

    const sq1 = '<svg width="20" height="20"> <rect width="20" height="20" style="fill:' + color1 +';stroke-width:3;stroke:rgb(0,0,0)" /> </svg>'

    const sq2 = '<svg width="20" height="20"> <rect width="20" height="20" style="fill:' + color2 +';stroke-width:3;stroke:rgb(0,0,0)" /> </svg>'
    console.log(color1);

    document.getElementById('square1').innerHTML = sq1;
    console.log(document.getElementById('square1'));
    document.getElementById('square2').innerHTML = sq2;
  });
