

function Car(size) {
    RigidBody.call(this);
    
    this.size = size;
    this.wheels = [];
    var halfSize = geoff.math.mulVectorScalar(this.size, 0.5);
    //front wheels
    this.wheels[0] = new geoff.Wheel([halfSize[0], halfSize[1]], 0.5);
    this.wheels[1] = new geoff.Wheel([-halfSize[0], halfSize[1]], 0.5);
    //rear wheels
    this.wheels[2] = new geoff.Wheel([halfSize[0], -halfSize[1]], 0.5);
    this.wheels[3] = new geoff.Wheel([-halfSize[0], -halfSize[1]], 0.5);
    
    this.frictionCoeff = 0.3 / 20.5448542;
    this.frontalArea = 1.9; // m^2
    this.airDensity = 1.29; // kg/m^3
}
inherit(Car, RigidBody);

