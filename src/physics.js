function RigidBody() {
    // linear properties
    this.position = vec2.create(0.0, 0.0);
    this.velocity = vec2.create(0.0, 0.0);
    this.force = vec2.create(0.0, 0.0);
    this.mass = 1.0;
    // angular properties
    this.angle = 0.0;
    this.angularVelocity = 0.0;
    this.torque = 0.0;
    this.angularInertia = 1.0;
};

Object.extend(RigidBody.prototype, {
    setAngle: function(angle) {
        const pi2 = Math.PI * 2;
        angle = angle % pi2;
        if (angle < 0.0) {
            angle += pi2;
        }
        this.angle = angle;
    },
    addForce: function(force, offset){
        this.force = vec2.add(this.force, force);
        this.torque += vec2.cross(offset, force);
    },
    update: function(timeStep) {
        var acceleration = vec2.scale(this.force, this.mass);
        this.velocity = vec2.add(this.velocity, vec2.scale(timeStep, acceleration));
        this.position = vec2.add(this.position, vec2.scale(timeStep, this.velocity));
        this.force = vec2.create(0.0, 0.0);
        
        var angularAccel = this.torque / this.angularInertia;
        this.angularVelocity += angularAccel * timeStep;
        this.setAngle(this.angle + this.angularVelocity * timeStep);
        this.torque = 0.0;
    },
    velocityOfPoint: function(worldOffset) {
        var tangent = vec2.create(-worldOffset[1], worldOffset[0]);
        return vec2.add(vec2.scale(tangent, this.angularVelocity), this.velocity);
    },
    relativeToWorld: function(relative) {
        return vec2.rotate(relative, this.angle);
    },
    worldToRelative: function(world) {
        return vec2.rotate(world, -this.angle);
    }
});
