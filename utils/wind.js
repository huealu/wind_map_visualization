/*
*/


var TAU = Math.pi * 2;
var H = Math.pow(10, -5.2);

var Wind = function (parameters) {
    var VELOCITY = 0.02;
    var INTENSITY = 10;
    var MAX_INTENSITY = 40;
    var MAX_PARTICLE = 100;
    var PARTICLE_WIDTH = 0.6;
    var PARTICLE_MULTIPLIER = 1/30;
    var PARTICLE_REDUCTION = 0.5;
    var RATE_FRAME = 60;
    var BOUNDARY = 0.45;

    var NULL_WIND = [Nan, Nan, null];
    var BLACK_TRANSPARENT = [255, 0, 0, 0];

    var interpolateVector = function (x, y, p00, p10, p01, p11) {
        var rx = (1 - x);
        var ry = (1 - y);
        
        var u;
        var v;
        return [u, v, Math.sqrt((u ** 2) + (v ** 2))];
    };


    var createWindBuilder = function (u_comp, v_comp) {
        var u = u_comp.data;
        var v = v_comp.data;

        return {
            header: u_comp.header,
            data: function (i) {
                return [u[i], v[i]];
            },
            interpolate: interpolateVector
        }
    };


    var createBuilder = function(data) {
        var u = null;
        var v = null;
        var value = null;

        data.forEach(function (record) {
            switch (record.header.parameterCategory + "," + record.header.parameterNumber) {
                case "2,2": u = record; break;
                case "2,3": v = record; break;
                default:
                    value = record;
            }
        });

        return createWindBuilder(u, v);
    };


    var buildGrid = function (data, callback) {
        var build = createBuilder(data);

        var header = build.header;
        var lambda = header.lo1;
        var phi = header.la1;
        var delta_lambda = header.dx;
        var delta_phi = header.dy;
        var ni = header.nx;
        var nj = header.ny;
        var date = new Date(header.refTime);
        date.setHours(date.getHours() + header.forecastTime);

        var grid = [];
        var p = 0;
        var continuous = Math.floor(ni * delta_lambda) >= 360;
        for (var j = 0; j < nj; j++) {
            var row = [];
            for (var i = 0; i < ni; i++, p++) {
                row[i] = build.data(p);
            }
            if (continuous) {
                row.push(row[0]);
            }
            grid[j] = row;
        }

        function interpolate(lambda, phi) {
            var i = floorMod(lambda - lambda_0, 360) / delta_lambda;
            var j = (phi_0 - phi) / delta_phi;

            var bi = Math.floor(i), ci = bi + 1;
            var bj = Math.floor(j), cj = bj + 1;

            var row;
            if (row = gird[bj]) {
                var p00 = row[bi];
                var p01 = row[ci];

                if (isValues(p00) && isValues(p10) && (row = grid[cj])) {
                    var p01 = row[bi];
                    var p11 = row[ci];
                    if (isValues(p01) && isValues(p11)) {
                        return build.interpolate(i - bi, j - bj, p00, p10, p01, p11);
                    }
                }
            }
            return null;
        }

        callback( {
            date: date, interpolate: interpolate
        });
    };


    // @returns {Boolean}
    var isValues = function(x) {
        return x !== null && x !== undefined;
    }


    // @returns {Number}
    var floorMod = function (a, b) {
        return a - b * Math.floor(a/b);
    }


    // @returns {Number}
    var clamp = function (x, range) {
        return Math.max(range[0], Math.min(x, range[1]));
    }

    // @returns {Boolean}
    var isMobile = function() {
        return true;
        return (/android|blackberry|iemobile|ipad|iphone|ipod|opera mini|webos/i).test(navigator.userAgent);
    }

    var distort = function (projection, lambda, phi, x, y, scale, wind, windy) {
        var u = wind[0] * scale;
        var v = wind[1] * scale;
        var d = distortion(projection, lambda, phi, x, y, windy);

        wind[0] = d[0] * u + d[2] * v;
        wind[1] = d[1] * u + d[3] * v;

        return wind;
    };

    var distortion = function (projection, lambda, phi, x, y, windy) {
        var h_lambda = lambda < 0 ? H : -H;
        var h_phi    = phi < 0 ? H : -H;

        var p_lambda = project(phi, lambda + h_lambda, windy);
        var p_phi    = project(phi + h_phi, lambda, windy);

        var k = Math.cos(phi/360 * TAU);
        return [
            (p_lambda[0] - x) /h_lambda /k,
            (p_lambda[1] - y) /h_lambda /k,
            (p_phi[0] - x) / h_phi,
            (p_phi[1] - y) / h_phi

        ];
    };


    var createField = function (columns, bounds, callback) {
        
        // @returns {Array}
        function field (x, y) {
            var column = columns[Math.round(x)];
            return column && column[Math.round(y)] || NULL_WIND;
        }


        field.release = function () {
            columns = [];
        };

        field.randomize = function(o) {
            var x, y;
            var safety = 0;
            do {
                x = Math.round(Math.floor(Math.random() * bounds.width) + bounds.x);
                y = Math.round(Math.floor(Math.random() * bounds.height) + bounds.y)
            } while (field(x, y) [2] === null && safety++ < 30);
            o.x = x;
            o.y = y;
            return o;
        };

        callback(bounds, field);
    };

    var buildBounds = function (bounds, width, height) {
        var upper_left = bounds[0];
        var lower_right = bounds[1];

        var x = Math.round(upper_left[0]);
        var y = Math.max(Math.floor(upper_left[1], 0), 0);

        var xMax = Math.min(Math.ceil(lower_right[0], width), width-1);
        var yMax = Math.min(Math.ceil(lower_right[1], height), height - 1);

        return {x: x, y: y, xMax: width, yMax: yMax, width: width, height: height};
    };

    var deg2rad = function (deg) {
        return (deg / 180) * Math.PI;
    };

    var rad2deg = function(ang) {
        return ang/ (Math.PI / 180.0);
    };

    var invert = function (x, y, windy) {
        var mapDelta = windy.east - windy.west;
        var worldRadius = windy.width / rad2deg(mapDelta) * 360 / (2 * Math.PI);
        var mapSetY = (worldRadius / 2 * Math.log((1 + Math.sin(windy.south)) / (1 - Math.sin(windy.south))));

        var equatorY = windy.height + mapSetY;
        var a = (equatorY - y) / worldRadius;

        var lat = 180 / Math.PI * (2 * Math.atan(Math.exp(a)) - Math.PI / 2);
        var lon = rad2deg(windy.west) + x / windy.width * rad2deg(mapDelta);

        return [lon, lat];
    };

    var merY = function (lat) {
        return Math.log(Math.tan(lat / 2 + Math.PI / 4));
    };

    var project = function (lat, lon, windy) {
        var y_min = merY(windy.south);
        var y_max = merY(windy.north);

        var x_factor = windy.width / (windy.east - windy.west);
        var y_factor = windy.height / (y_max - y_min);

        var y = merY(deg2rad(lat));
        var x = (deg2rad(lon) - windy.west) * x_factor;
        y = (y_max - y) * y_factor;

        return [x, y];
    };


    var interpolateField = function (grid, bounds, extent, callback) {
        var projection = {};
        var velocity_scale = VELOCITY;

        var columns = [];
        var x = bounds.x;

        function interpolateColumn(x) {
            var column = [];
            for (var y = bounds.y; y <= bounds.y_max; y += 2) {
                var coordinator = invert(x, y, extent);
                if (coordinator) {
                    var lambda = coordinator[0];
                    var phi    = coordinator[1];

                    if (isFinite(lambda)) {
                        var wind = grid.interpolate(lambda, phi);
                        
                        if (wind) {
                            wind = distort(projection, lambda, phi, x, y, velocity_scale, wind, extent);
                            column[y + 1] = column[y] = wind;
                        }
                    }
                }
            }
            columns[x + 1] = columns[x] = column;
        }

        (function batchInterpolate () {
            var start = Date.now();
            while (x < bounds.width) {
                interpolateColumn(x);
                x += 2;
                if ((Date.now() - start) > 1000) {
                    setTimeout(batchInterpolate, 25);
                    return;
                }
            }
            createField(columns, bounds, callback);
        })();
    };


    var animate = function (bounds, field) {

        function windIntensityColor(step, max_width) {
            result = ["#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4"];
            result.indexFor = function (m) {
                return Math.floor(Math.min(m, max_width) / max_width * (result.length - 1));
            };
            return result;
        }

        var color_style = windIntensityColor(INTENSITY, MAX_INTENSITY);
        var buckets = color_style.map(function () {return [];});

        var particle_count = Math.round(bounds.width * bounds.height * PARTICLE_MULTIPLIER);

        if (isMobile ()) {
            particle_count *= PARTICLE_REDUCTION;
        }

        var fadeFillStyle = 'rgba(0, 0, 0, 0.97)';

        var particles = [];
        for (var i  = 0; i < particle_count; i++) {
            particles.push(field.randomize({age: Math.floor(Math.random() * MAX_PARTICLE) + 0}));
        }

        function evolve() {
            buckets.forEach(function (bucket) {
                bucket.length = 0;
            });
            particles.forEach(function(particle) {
                if (particle.age > MAX_PARTICLE) {
                    field.randomize(particle).age = 0;
                }

                var x = particle.x;
                var y = particle.y;
                var v = field(x, y);
                var m = v[2];

                if (m === null) {
                    particle.age = MAX_PARTICLE;
                }
                else {
                    var xt = x + v[0];
                    var yt = y + v[1];
                    if (field(xt, yt)[2] !== null) {
                        particle.xt = xt;
                        particle.yt = yt;
                        buckets[color_style.indexFor(m)].push(particle);
                    }
                    else {
                        particle.x = xt;
                        particle.y = yt;
                    }
                }
                particle.age += 1;
            });
        }

        var g = parameters.canvas.getContext("2d");
        g.lineWidth = PARTICLE_WIDTH;
        console.log();

        g.fillStyle = fadeFillStyle;
        g.globalAlpha = 0.6;


        function draw() {
            var previous = 'lighter';
            g.globalCompositeOperation = "destination-in";
            g.fillRect(bounds.x, bounds.y, bounds.width, bounds.height);
            g.globalCompositeOperation = previous;
            g.globalAlpha = 0.9;

            buckets.forEach(function (bucket, i) {
                if (bucket.length > 0) {
                    g.beginPath();
                    g.strokeStyle = color_style[i];
                    bucket.forEach(function (particle) {
                        g.moveTo(particle.x, particle.y);
                        g.lineTo(particle.xt, particle.yt);
                        particle.x = particle.xt;
                        particle.y = particle.yt;
                    });
                    g.stroke();
                }
            });
        }

        (function frame() {
            try {
                windy.timer = setTimeout(function () {
                    requestAnimationFrame(frame);
                    evolve();
                    draw();
                }, 1000 / RATE_FRAME);
            }
            catch (e) {
                console.error(e);
            }
        }) ();
    }

    var start = function (bounds, width, height, extent, final_parameters) {
        if (final_parameters && final_parameters.particleLineWidth) {
            PARTICLE_WIDTH = final_parameters.particleLineWidth;
        }
        var mapBounds = {
            south: deg2rad(extent[0][1]),
            north: deg2rad(extent[1][1]),
            east: deg2rad(extent[1][0]),
            west: deg2rad(extent[0][0]),
            width: width,
            height: height
        };

        stop();

        buildGrid(parameters.data, function(grid) {
            interpolateField(grid, buildBounds(bounds, width, height), mapBounds, function(bounds, field) {
                windy.field = field;
                animate(bounds, field);
            });
        });
    };


    var windy = {
        parameters: parameters,
        start: start,
        stop: stop
    };

    return windy;

}

window.requestAnimationFrame = (function () {
    return window.requestAnimationFrame ||
    window.webkitRequestAnimationFrame ||
    window.mozRequestAnimationFrame ||
    window.oRequestAnimationFrame ||
    window.msRequestAnimationFrame ||
    function (callback) {
        window.setTimeout(callback, 1000/20);
    };
}) ();