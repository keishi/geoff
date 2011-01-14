function inherit(childClass, parentClass) {
    var tempClass = function() {};
    tempClass.prototype = parentClass.prototype;
    childClass.prototype = new tempClass();
    childClass.prototype.constructor = childClass;
};