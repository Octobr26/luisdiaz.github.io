/*globals paper, console, $ */
/*jslint nomen: true, undef: true, sloppy: true */

/*

@licstart  The following is the entire license notice for the
JavaScript code in this page.

Copyright (C) 2015 david ha, otoro.net, otoro labs

The JavaScript code in this page is free software: you can
redistribute it and/or modify it under the terms of the GNU
General Public License (GNU GPL) as published by the Free Software
Foundation, either version 3 of the License, or (at your option)
any later version.  The code is distributed WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU GPL for more details.

As additional permission under GNU GPL version 3 section 7, you
may distribute non-source (e.g., minimized or compacted) forms of
that code without the copy of the GNU GPL normally required by
section 4, provided you include this license notice and a URL
through which recipients can access the Corresponding Source.


@licend  The above is the entire license notice
for the JavaScript code in this page.
*/

var CREATURE = CREATURE || {};

///////////////////////////
// start main code here //
///////////////////////////
CREATURE.Main = function (mobileMode) {
    // put paper.js in local environment, and setup canvas and tool (for events)
    this.paper = new paper.PaperScope();
    var canvasName = "topCanvas";
    this.paper.setup(canvasName);

    with (this.paper) {
        var DETAILED_MODE = false; // true but uses more computer time.

        if (mobileMode) {
            console.log("in mobile mode, switching off detail");
            DETAILED_MODE = false;
        } else {
            console.log("in desktop mode, do nothing now");
        }

        var tool = new Tool();

        // hard coded settings for different types of creatures:
        // (in the future, dynamically change the size of creatures depending on resolution)
        var foodCreatureSettings = {
            creatureType: "food",
            radius: 16,
            maxSpeed: 2,
            maxForce: 0.05,
            nMembranePoints: 8,
            bodyColor: "#5AC74E",
            bounciness: 0.11,
            smoothBody: true
        };
        var hunterCreatureSettings = {
            creatureType: "hunter",
            radius: 48,
            maxSpeed: 5,
            maxForce: 0.15,
            nMembranePoints: 16,
            bodyColor: "#C83232",
            bounciness: 0.05,
            smoothBody: true
        };

        // creature object
        var Creature = function (position, creatureSettings) {
            var i = 0,
                thetaStep = 0;
            creatureSettings = creatureSettings || foodCreatureSettings;

            this.target = Point.random().multiply(view.size);
            this.position = position.clone();
            this.velocity = Point.random().multiply(2.0).add([-1.0, -1.0]);
            this.acceleration = new Point(0, 0);
            this.maxSpeed = creatureSettings.maxSpeed * getRandom(0.75, 1.25);
            this.maxForce = creatureSettings.maxForce;
            this.radius = creatureSettings.radius;
            this.bounciness = creatureSettings.bounciness * getRandom(0.8, 1.2);
            this.nMembranePoints = creatureSettings.nMembranePoints;
            this.cid = getNewCreatureID();
            this.creatureType = creatureSettings.creatureType;
            this.isAlive = true; // if set to dead, will remove from creatures stack next frame.

            this.membraneDeviation = new Array(this.nMembranePoints);
            this.membraneLocX = new Array(this.nMembranePoints);
            this.membraneLocY = new Array(this.nMembranePoints);
            thetaStep = 360 / this.nMembranePoints;
            for (i = 0; i < this.nMembranePoints; i += 1) {
                this.membraneDeviation[i] = 0;
                this.membraneLocX[i] = fastSin(i * thetaStep);
                this.membraneLocY[i] = fastCos(i * thetaStep);
            }

            this.bodyColor = creatureSettings.bodyColor;
            this.smoothBody = creatureSettings.smoothBody;
            this.createItems();
        };

        Creature.prototype.run = function () {
            this.separate();
            this.seek();
            this.update();
            this.moveHead();
        };

        Creature.prototype.createItems = function () {
            var i = 0,
                r = this.radius;

            this.eye = [];
            this.pupil = [];

            for (i = 0; i < 2; i++) {
                if (this.radius > 8) {
                    // only draw eye pupil if radius is large enough
                    this.pupil[i] = new Path.Ellipse({
                        center: [0, 0],
                        radius: r / 6,
                        fillColor: "black"
                    }).sendToBack();
                }

                this.eye[i] = new Path.Ellipse({
                    center: [0, 0],
                    radius: r / 3,
                    fillColor: "white",
                    strokeColor: "black"
                }).sendToBack();
            }

            this.head = new Path();

            this.head.strokeColor = "black";
            i = this.nMembranePoints - 1;
            for (; i >= 0; i--) {
                this.head.add(
                    new Point(
                        r * this.membraneLocX[i],
                        r * this.membraneLocY[i]
                    )
                );
            }
            this.head.closed = true;
            this.head.fillColor = this.bodyColor;

            if (DETAILED_MODE) {
                this.head.opacity = getRandom(0.4, 0.9);
                //this.head.fillColor.hue = getRandom(0.5, 1.0);
            }

            this.head.sendToBack();
        };
        Creature.prototype.removeItems = function () {
            this.eye[0].remove();
            this.eye[1].remove();
            if (this.pupil[0]) {
                this.pupil[0].remove();
                this.pupil[1].remove();
            }
            this.head.remove();
        };

        Creature.prototype.animateItems = function () {
            var i = this.nMembranePoints - 1,
                deviation = 0.0,
                r = this.radius,
                eyeIndex,
                eyeDist = 0.6,
                len,
                antilen,
                locX,
                locY;
            for (; i >= 0; i--) {
                deviation = this.membraneDeviation[i];
                deviation += r * this.bounciness * getRandom(-1.0, 1.0);
                deviation *= 0.95;
                this.membraneDeviation[i] = deviation;
                len = r + deviation;
                locX = this.membraneLocX[i];
                locY = this.membraneLocY[i];
                this.head.segments[i].point.x = len * locX;
                this.head.segments[i].point.y = len * locY;
            }
            if (SMOOTH_MODE && this.smoothBody) {
                this.head.smooth();
            }

            for (i = 0; i < 2; i++) {
                eyeIndex = Math.round(
                    this.nMembranePoints * (1 - (0.625 + 0.25 * i))
                );
                len = r + this.membraneDeviation[eyeIndex] * 1;
                antilen = r + this.membraneDeviation[eyeIndex + 1 - 2 * i] * 1;
                locX = this.membraneLocX[eyeIndex];
                locY = this.membraneLocY[eyeIndex];
                this.eye[i].position = [
                    (len * locX * eyeDist * len) / r,
                    (len * locY * eyeDist * antilen) / r
                ];
                this.eye[i].radius = [len / 3, len / 3];
                if (this.pupil[i]) {
                    // if creature has eye pupil
                    this.pupil[i].position = [
                        (len * locX * eyeDist * len) / r,
                        (len * locY * eyeDist * antilen) / r
                    ];
                    this.pupil[i].radius = [len / 6, len / 6];
                }
            }
        };

        Creature.prototype.moveHead = function () {
            var i;
            this.animateItems();
            this.head.position = this.position;
            this.head.rotate(this.velocity.angle, this.position);
            for (i = 0; i < 2; i++) {
                this.eye[i].position.x += this.position.x;
                this.eye[i].position.y += this.position.y;
                this.eye[i].rotate(this.velocity.angle, this.position);
                if (this.pupil[i]) {
                    // if creature has eye pupil
                    this.pupil[i].position.x += this.position.x;
                    this.pupil[i].position.y += this.position.y;
                    this.pupil[i].rotate(this.velocity.angle, this.position);
                }
            }
            //this.targetMarker.position = this.target;
        };

        // A method that calculates a steering vector towards a target
        Creature.prototype.steer = function () {
            var steer,
                desired = this.target.subtract(this.position);
            var distance = desired.length;

            if (distance < this.radius) {
                // slowly go towards target (removed)
                desired.length = this.maxSpeed * (distance / this.radius);
            } else {
                desired.length = this.maxSpeed;
            }
            //desired.length = this.maxSpeed;
            steer = desired.subtract(this.velocity);
            steer.length = Math.min(this.maxForce, steer.length);

            return steer;
        };

        Creature.prototype.seek = function () {
            var direction = this.getNearestEnemyDirection();
            var nudge;
            if (direction) {
                // set target to/away from closest enemy:
                this.target = this.position.add(direction);
            } else {
                // if no more food is available, then randomy wander:
                this.setTargetFromAction(getRandomInt(0, nActionState));
                // also nudge towards the center:
                nudge = view.center.subtract(this.target);
                nudge.length = this.radius / 1.5;
                this.target = this.target.add(nudge);
            }

            var separation = this.separate().multiply(3.0);
            this.acceleration = this.acceleration.add(separation);
            this.acceleration = this.acceleration.add(this.steer());
        };

        Creature.prototype.getNearestEnemyDirection = function () {
            var bestDirection = null;
            var bestDistance = LARGENUM;
            var i = creatures.length - 1;
            var other;
            var direction = null;
            var distance;
            for (; i >= 0; i--) {
                other = creatures[i];
                if (other.creatureType != this.creatureType && other.isAlive) {
                    direction = other.position.subtract(this.position);
                    distance = direction.length;
                    if (distance < bestDistance) {
                        bestDistance = distance;
                        bestDirection = direction;
                    }
                }
            }
            if (this.creatureType === "food") {
                bestDirection = bestDirection.multiply(-1.0); // move away from hunter if food.
            }
            return bestDirection;
        };

        // based off a number from 0 to nActionState (exclisive) return the point on the circle that will be the target for agent.
        Creature.prototype.setTargetFromAction = function (/* int */ action) {
            var angle = (360.0 * action) / nActionState;
            var direction = this.velocity.clone();
            var theTarget = new Point(
                fastCos(angle) * this.radius,
                fastSin(angle) * this.radius
            );
            direction.length = this.radius;
            this.target = theTarget.add(direction).add(this.position);
        };

        Creature.prototype.separate = function () {
            var desiredSeperation;
            var steer = new Point(0, 0);
            var count = 0;
            var i = creatures.length - 1;
            var other;
            var distVector;
            var distance;
            // For every creature in the system, check if it's too close
            for (; i >= 0; i--) {
                other = creatures[i];
                desiredSeperation = this.radius + other.radius;
                if (other.cid && other.cid != this.cid) {
                    distVector = this.position.subtract(other.position);
                    distance = distVector.length;
                    if (
                        distance <= desiredSeperation &&
                        distance > 0 &&
                        other.creatureType === this.creatureType
                    ) {
                        // Calculate vector pointing away from neighbor
                        steer = steer.add(distVector.normalize(1 / distance));
                        count++;
                    } else if (
                        this.creatureType === "hunter" &&
                        other.isAlive &&
                        other.creatureType != this.creatureType &&
                        distance <= Math.max(this.radius, other.radius)
                    ) {
                        // other creature died
                        other.isAlive = false;
                    }
                }
            }
            // Average -- divide by how many
            if (count > 0) {
                steer = steer.divide(count);
                // Implement Reynolds: Steering = Desired - Velocity
                steer.length = this.maxSpeed;
                steer = steer.subtract(this.velocity);
                steer.length = Math.min(steer.length, this.maxForce);
            }
            return steer;
        };

        Creature.prototype.update = function () {
            // Update velocity
            this.velocity = this.velocity.add(this.acceleration);
            // Limit speed (vector#limit?)
            this.velocity.length = Math.min(
                this.maxSpeed,
                this.velocity.length
            );
            this.position = this.position.add(this.velocity);

            //this.position.x = (this.position.x+view.size.width) % view.size.width;
            //this.position.y = (this.position.y+view.size.height) % view.size.height;

            // Reset acceleration to 0 each cycle
            this.acceleration.x = 0;
            this.acceleration.y = 0;
        };

        tool.onMouseDown = function (event) {
            addFoodAtMousePoint(event.point);
        };

        tool.onMouseDrag = function (event) {
            addFoodAtMousePoint(event.point);
        };

        var addFoodAtMousePoint = function (point) {
            var creature = getCreatureAtPosition(point);
            if (!creature) {
                addFoodCreature(point);
            } else {
                // only drag big hunter creatures
                if (creature.creatureType === "hunter") {
                    creature.position = point;
                }
            }
        };

        var addFoodCreature = function (point) {
            if (creatures.length < MAX_CREATURES) {
                creatures.push(new Creature(point, foodCreatureSettings));
            }
        };

        var draw = function (event) {
            var Nstart = 480,
                Nd = 240,
                percent;

            if (event.count >= Nstart && event.count <= Nstart + Nd) {
                percent = (event.count - Nstart) / Nd;

                title.opacity = 0.5 * percent;

                title.style.fontSize =
                    titleFontSize *
                    (percent +
                        percent *
                            fastSin(Math.exp(percent * percent) * 4800) *
                            Math.exp(-4.0 * percent));
            }

            if (
                event.count > Nstart + Nd * 2 &&
                event.count <= Nstart + Nd * 2.5
            ) {
                title.opacity *= 0.95;
            }

            if (
                event.count >= Nstart + Nd * 2.5 &&
                event.count <= Nstart + Nd * 5
            ) {
                title.opacity = Math.min(0.5, title.opacity + 0.005);
                title.content = "scroll->down();";
            }

            if (event.count > Nstart + Nd * 5) {
                title.opacity *= 0.95;
            }

            if (event.count === 60) {
                // if the average time of first 120 frames is < 40fps
                //switch off detailed mode
                if (event.time > 3.5) {
                    console.log("slow browser mode started");
                    DETAILED_MODE = false;
                    setAllCreatureOpacity(1.0);
                    MAX_CREATURES = DETAILED_MODE ? 14 : 8;
                    //SMOOTH_MODE = false;
                } else {
                    console.log("fast browser confirmed");
                    SMOOTH_MODE = true;
                }
            }

            // update fps:
            if ((event.count + 1) % 60 === 0) {
                /*
			fps_data.content = Math.round(10 / event.delta) / 10;
			*/
                fps_data.content = "";
            }

            var len = creatures.length,
                i = len - 1;

            for (; i >= 0; i--) {
                if (creatures[i].isAlive) {
                    // make each creature run according to their logic:
                    creatures[i].run();
                } else {
                    // remove dead creatures
                    creatures[i].removeItems();
                    creatures.splice(i, 1);
                }
            }
        };

        // useful helper functions:

        var fastSin = function (xDeg) {
            var deg = Math.round(xDeg);
            if (deg >= 0) {
                return sinTable[deg % 360];
            }
            return -sinTable[-deg % 360];
        };

        var fastCos = function (xDeg) {
            var deg = Math.round(Math.abs(xDeg));
            return cosTable[deg % 360];
        };

        var getRandom = function (min, max) {
            return Math.random() * (max - min) + min;
        };

        var getRandomInt = function (min, max) {
            return Math.floor(Math.random() * (max - min)) + min;
        };

        var getNewCreatureID = function () {
            creature_id += 1;
            return creature_id;
        };

        var getCreatureAtPosition = function (position) {
            var result = null;
            var i = creatures.length - 1;
            var distVector;
            var distance;
            // For every creature in the system, check if it's too close
            for (; i >= 0; i--) {
                distVector = position.subtract(creatures[i].position);
                distance = distVector.length;
                if (distance <= creatures[i].radius) {
                    result = creatures[i];
                    break;
                }
            }
            return result;
        };

        var setAllCreatureOpacity = function (opacityLevel) {
            var i = creatures.length - 1;
            for (; i >= 0; i--) {
                creatures[i].head.opacity = opacityLevel;
            }
        };

        // INIT:

        var creatures = [];
        var i = 0;
        var cosTable = new Array(360);
        var sinTable = new Array(360);
        var PI = Math.PI;
        var nActionState = 360;
        var position;
        var creature_id = 0; // create unique identifier
        var LARGENUM = 1000000000000000;
        var screenRadius = view.center.length;
        var randomAngle;
        var MAX_CREATURES = DETAILED_MODE ? 16 : 12;
        var SMOOTH_MODE = true;

        // pre compute sine and cosine values to the nearest degree
        for (i = 0; i < 360; i++) {
            cosTable[i] = Math.cos((i / 360) * 2 * PI);
            sinTable[i] = Math.sin((i / 360) * 2 * PI);
        }

        // Add the creatures:
        for (i = 0; i < 3; i++) {
            randomAngle = getRandom(0, 360);
            position = new Point(
                fastCos(randomAngle) * screenRadius * 2,
                fastSin(randomAngle) * screenRadius * 2
            );
            creatures.push(
                new Creature(position.add(view.center), hunterCreatureSettings)
            );
        }

        for (i = MAX_CREATURES - creatures.length - 1 - 4; i >= 0; i--) {
            position = Point.random().multiply(view.size);
            creatures.push(new Creature(view.center, foodCreatureSettings));
        }

        // title:

        var title = new PointText(view.size);
        var titleFontSize = screenRadius / 10;
        title.position.y -= 32;
        title.position.x -= 32;
        title.content = "바나나 우유"; //'{ banana-uyu.net };';
        title.style = {
            fontFamily: "Inter",
            fontWeight: "normal",
            fontSize: 0,
            //fillColor: '#1664ff',
            fillColor: "#1664ff",
            justification: "right"
        };

        // fps handler
        var fps_data = new PointText(8, 24);
        fps_data.content = "";
        fps_data.style = {
            fontFamily: "Courier New",
            fontWeight: "normal",
            fontSize: 16,
            //fillColor: '#1664ff',
            fillColor: "#cacaca",
            justification: "left"
        };
        if (DETAILED_MODE) {
            fps_data.opacity = 1.0;
        }

        window.onload = function () {
            view.onFrame = function (event) {
                if (drawCreatureStatus === true) {
                    draw(event);
                }
            };

            view.onResize = function (event) {
                // Whenever the window is resized, recenter the path:
                console.log("resized window.");
                title.point = view.size.clone();
                title.position.y -= 32;
                title.position.x -= 32;
                title.style.justification = "right";
            };
        };
    }
};

//////////////////
// end ///////////
//////////////////
