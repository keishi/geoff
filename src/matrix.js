var matrixArrayType = Array;
if(typeof Float32Array != 'undefined') {
	glMatrixArrayType = Float32Array;
}

var vec2 = {
    create: function(x, y) {
        var dest = new matrixArrayType(2);
		if (typeof x != "undefined") {
            dest[0] = x;
    		dest[1] = y;
        }
        return dest;
    },
    copy: function(v){
        var dest = new matrixArrayType(2);
        dest[0] = v[0];
    	dest[1] = v[1];
    	return dest;
    },
    add: function(a, b, dest){
        dest = dest || new matrixArrayType(2);
        dest[0] = a[0] + b[0];
    	dest[1] = a[1] + b[1];
        return dest;
    },
    subtract: function(a, b, dest){
        dest = dest || new matrixArrayType(2);
        dest[0] = a[0] - b[0];
    	dest[1] = a[1] - b[1];
        return dest;
    },
    negate: function(v, dest){
        dest = dest || new matrixArrayType(2);
        dest[0] = -v[0];
    	dest[1] = -v[1];
        return dest;
    },
    scale: function(v, f, dest){
        dest = dest || new matrixArrayType(2);
        dest[0] = f * v[0];
    	dest[1] = f * v[1];
        return dest;
    },
    normalize: function(v, dest){
        dest = dest || new matrixArrayType(2);
        var x = v[0], y = v[1];
    	var f = 1 / Math.sqrt(x*x + y*y);
    	dest[0] = x * f;
    	dest[1] = y * f;
    	return dest;
    },
    cross: function(a, b){
    	return a[0] * b[1] - a[1] * b[0];
    },
    length: function(v){
        var x = v[0], y = v[1];
        return Math.sqrt(x*x + y*y);
    },
    dot: function(a, b){
        return a[0] * b[0] + a[1] * b[1];
    },
    direction: function(attribute){
        dest = dest || new matrixArrayType(2);
        var x = a[0] - b[0];
        var y = a[1] - b[1];
    	var f = 1 / Math.sqrt(x*x + y*y);
    	dest[0] = x * f;
    	dest[1] = y * f;
    	return dest;
    },
    rotate: function(v, angle, dest){
        dest = dest || new matrixArrayType(2);
        var c = Math.cos(angle);
        var s = Math.sin(angle);
        dest[0] = v[0] * c - v[1] * s;
        dest[1] = v[0] * s + v[1] * c;
        return dest;
    }
};

var vec3 = {
    create: function(x, y, z) {
        var dest = new matrixArrayType(3);
        if (typeof x != "undefined") {
            dest[0] = x;
    		dest[1] = y;
    		dest[2] = z;
        }
        return dest;
    },
    copy: function(v){
        var dest = new matrixArrayType(3);
        dest[0] = v[0];
    	dest[1] = v[1];
    	dest[2] = v[2];
    	return dest;
    },
    add: function(a, b, dest){
        dest = dest || new matrixArrayType(3);
        dest[0] = a[0] + b[0];
    	dest[1] = a[1] + b[1];
    	dest[2] = a[2] + b[2];
        return dest;
    },
    subtract: function(a, b, dest){
        dest = dest || new matrixArrayType(3);
        dest[0] = a[0] - b[0];
    	dest[1] = a[1] - b[1];
    	dest[2] = a[2] - b[2];
        return dest;
    },
    negate: function(v, dest){
        dest = dest || new matrixArrayType(3);
        dest[0] = -v[0];
    	dest[1] = -v[1];
    	dest[2] = -v[2];
        return dest;
    },
    scale: function(v, f, dest){
        dest = dest || new matrixArrayType(3);
        dest[0] = f * v[0];
    	dest[1] = f * v[1];
    	dest[2] = f * v[2];
        return dest;
    },
    normalize: function(v, dest){
        dest = dest || new matrixArrayType(3);
        var x = v[0], y = v[1], z = v[2];
    	var f = 1 / Math.sqrt(x*x + y*y + z*z);
    	dest[0] = x * f;
    	dest[1] = y * f;
    	dest[2] = z * f;
    	return dest;
    },
    cross: function(a, b, dest){
        dest = dest || new matrixArrayType(3);
        var ax = a[0], ay = a[1], az = a[2];
        var bx = b[0], by = b[1], bz = b[2];
    	dest[0] = ay * bz - az * by;
    	dest[1] = az * bx - ax * bz;
    	dest[2] = ax * by - ay * bx;
    	return dest;
    },
    length: function(v){
        var x = v[0], y = v[1], z = v[2];
        return Math.sqrt(x*x + y*y + z*z);
    },
    dot: function(a, b){
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    },
    direction: function(attribute){
        dest = dest || new matrixArrayType(3);
        var x = a[0] - b[0];
        var y = a[1] - b[1];
        var z = a[2] - b[2];
    	var f = 1 / Math.sqrt(x*x + y*y + z*z);
    	dest[0] = x * f;
    	dest[1] = y * f;
    	dest[2] = z * f;
    	return dest;
    }
};

var mat3 = {
    create: function(m11, m12, m13, m21, m22, m23, m31, m32, m33) {
        dest = dest || new matrixArrayType(9);
    	if (typeof m11 != "undefined") {
    		dest[0] = m11;
    		dest[1] = m12;
    		dest[2] = m13;
    		dest[3] = m21;
    		dest[4] = m22;
    		dest[5] = m23;
    		dest[6] = m31;
    		dest[7] = m32;
    		dest[8] = m33;
    	}
    	return dest;
    },
    copy: function(m) {
        dest = dest || new matrixArrayType(9);
		dest[0] = m[0];
		dest[1] = m[1];
		dest[2] = m[2];
		dest[3] = m[3];
		dest[4] = m[4];
		dest[5] = m[5];
		dest[6] = m[6];
		dest[7] = m[7];
		dest[8] = m[8];
    	return dest;
    },
    identity: function(dest){
        dest = dest || new matrixArrayType(9);
        dest[0] = 1;
    	dest[1] = 0;
    	dest[2] = 0;
    	dest[3] = 0;
    	dest[4] = 1;
    	dest[5] = 0;
    	dest[6] = 0;
    	dest[7] = 0;
    	dest[8] = 1;
    	return dest;
    },
    toMat4: function(mat, dest){
        dest = dest || new matrixArrayType(9);
        dest[0] = mat[0];
    	dest[1] = mat[1];
    	dest[2] = mat[2];
    	dest[3] = 0;

    	dest[4] = mat[3];
    	dest[5] = mat[4];
    	dest[6] = mat[5];
    	dest[7] = 0;

    	dest[8] = mat[6];
    	dest[9] = mat[7];
    	dest[10] = mat[8];
    	dest[11] = 0;

    	dest[12] = 0;
    	dest[13] = 0;
    	dest[14] = 0;
    	dest[15] = 1;

    	return dest;
    }
};

var mat4 = {
    create: function(m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44) {
        dest = dest || new matrixArrayType(15);
    	if (typeof m11 != "undefined") {
    		dest[0] = m11;
    		dest[1] = m12;
    		dest[2] = m13;
    		dest[3] = m14;
    		dest[4] = m21;
    		dest[5] = m22;
    		dest[6] = m23;
    		dest[7] = m24;
    		dest[8] = m31;
    		dest[9] = m32;
    		dest[10] = m33;
    		dest[11] = m34;
    		dest[12] = m41;
    		dest[13] = m42;
    		dest[14] = m43;
    		dest[15] = m44;
    	}
    	return dest;
    },
    identity: function(dest){
        dest = dest || new matrixArrayType(15);
        dest[0] = 1;
    	dest[1] = 0;
    	dest[2] = 0;
    	dest[3] = 0;
    	dest[4] = 0;
    	dest[5] = 1;
    	dest[6] = 0;
    	dest[7] = 0;
    	dest[8] = 0;
    	dest[9] = 0;
    	dest[10] = 1;
    	dest[11] = 0;
    	dest[12] = 0;
    	dest[13] = 0;
    	dest[14] = 0;
    	dest[15] = 1;
    	return dest;
    },
    transpose: function(m, dest){
        dest = dest || new matrixArrayType(15);
        dest[0] = m[0];
        dest[1] = m[4];
        dest[2] = m[8];
        dest[3] = m[12];
        dest[4] = m[1];
        dest[5] = m[5];
        dest[6] = m[9];
        dest[7] = m[13];
        dest[8] = m[2];
        dest[9] = m[6];
        dest[10] = m[10];
        dest[11] = m[14];
        dest[12] = m[3];
        dest[13] = m[7];
        dest[14] = m[11];
        dest[15] = m[15];
        return dest;
    },
    determinant: function(mat){
        var m00 = mat[0], m01 = mat[1], m02 = mat[2], m03 = mat[3];
    	var m10 = mat[4], m11 = mat[5], m12 = mat[6], m13 = mat[7];
    	var m20 = mat[8], m21 = mat[9], m22 = mat[10], m23 = mat[11];
    	var m30 = mat[12], m31 = mat[13], m32 = mat[14], m33 = mat[15];

    	return	m30*m21*m12*m03 - m20*m31*m12*m03 - m30*m11*m22*m03 + m10*m31*m22*m03 +
    			m20*m11*m32*m03 - m10*m21*m32*m03 - m30*m21*m02*m13 + m20*m31*m02*m13 +
    			m30*m01*m22*m13 - m00*m31*m22*m13 - m20*m01*m32*m13 + m00*m21*m32*m13 +
    			m30*m11*m02*m23 - m10*m31*m02*m23 - m30*m01*m12*m23 + m00*m31*m12*m23 +
    			m10*m01*m32*m23 - m00*m11*m32*m23 - m20*m11*m02*m33 + m10*m21*m02*m33 +
    			m20*m01*m12*m33 - m00*m21*m12*m33 - m10*m01*m22*m33 + m00*m11*m22*m33;
    },
    inverse: function(mat, dest) {
        dest = dest || new matrixArrayType(15);
        
        var a00 = mat[0], a01 = mat[1], a02 = mat[2], a03 = mat[3];
    	var a10 = mat[4], a11 = mat[5], a12 = mat[6], a13 = mat[7];
    	var a20 = mat[8], a21 = mat[9], a22 = mat[10], a23 = mat[11];
    	var a30 = mat[12], a31 = mat[13], a32 = mat[14], a33 = mat[15];

    	var b00 = a00*a11 - a01*a10;
    	var b01 = a00*a12 - a02*a10;
    	var b02 = a00*a13 - a03*a10;
    	var b03 = a01*a12 - a02*a11;
    	var b04 = a01*a13 - a03*a11;
    	var b05 = a02*a13 - a03*a12;
    	var b06 = a20*a31 - a21*a30;
    	var b07 = a20*a32 - a22*a30;
    	var b08 = a20*a33 - a23*a30;
    	var b09 = a21*a32 - a22*a31;
    	var b10 = a21*a33 - a23*a31;
    	var b11 = a22*a33 - a23*a32;
    	
    	var invDet = 1/(b00*b11 - b01*b10 + b02*b09 + b03*b08 - b04*b07 + b05*b06);

    	dest[0] = (a11*b11 - a12*b10 + a13*b09) * invDet;
    	dest[1] = (-a01*b11 + a02*b10 - a03*b09) * invDet;
    	dest[2] = (a31*b05 - a32*b04 + a33*b03) * invDet;
    	dest[3] = (-a21*b05 + a22*b04 - a23*b03) * invDet;
    	dest[4] = (-a10*b11 + a12*b08 - a13*b07) * invDet;
    	dest[5] = (a00*b11 - a02*b08 + a03*b07) * invDet;
    	dest[6] = (-a30*b05 + a32*b02 - a33*b01) * invDet;
    	dest[7] = (a20*b05 - a22*b02 + a23*b01) * invDet;
    	dest[8] = (a10*b10 - a11*b08 + a13*b06) * invDet;
    	dest[9] = (-a00*b10 + a01*b08 - a03*b06) * invDet;
    	dest[10] = (a30*b04 - a31*b02 + a33*b00) * invDet;
    	dest[11] = (-a20*b04 + a21*b02 - a23*b00) * invDet;
    	dest[12] = (-a10*b09 + a11*b07 - a12*b06) * invDet;
    	dest[13] = (a00*b09 - a01*b07 + a02*b06) * invDet;
    	dest[14] = (-a30*b03 + a31*b01 - a32*b00) * invDet;
    	dest[15] = (a20*b03 - a21*b01 + a22*b00) * invDet;
    	
    	return dest;
    },
    frustum: function(left, right, bottom, top, near, far, dest){
        dest = dest || new matrixArrayType(15);
        var rl = (right - left);
    	var tb = (top - bottom);
    	var fn = (far - near);
    	dest[0] = (near*2) / rl;
    	dest[1] = 0;
    	dest[2] = 0;
    	dest[3] = 0;
    	dest[4] = 0;
    	dest[5] = (near*2) / tb;
    	dest[6] = 0;
    	dest[7] = 0;
    	dest[8] = (right + left) / rl;
    	dest[9] = (top + bottom) / tb;
    	dest[10] = -(far + near) / fn;
    	dest[11] = -1;
    	dest[12] = 0;
    	dest[13] = 0;
    	dest[14] = -(far*near*2) / fn;
    	dest[15] = 0;
    	return dest;
    },
    perspective: function(fovy, aspect, near, far, dest){
        var top = near*Math.tan(fovy*Math.PI / 360.0);
    	var right = top*aspect;
    	return mat4.frustum(-right, right, -top, top, near, far, dest);
    },
    ortho: function(left, right, bottom, top, near, far, dest){
        dest = dest || new matrixArrayType(15);
        var rl = (right - left);
    	var tb = (top - bottom);
    	var fn = (far - near);
    	dest[0] = 2 / rl;
    	dest[1] = 0;
    	dest[2] = 0;
    	dest[3] = 0;
    	dest[4] = 0;
    	dest[5] = 2 / tb;
    	dest[6] = 0;
    	dest[7] = 0;
    	dest[8] = 0;
    	dest[9] = 0;
    	dest[10] = -2 / fn;
    	dest[11] = 0;
    	dest[12] = -(left + right) / rl;
    	dest[13] = -(top + bottom) / tb;
    	dest[14] = -(far + near) / fn;
    	dest[15] = 1;
    	return dest;
    },
    ortho: function(left, right, bottom, top, near, far, dest){
        dest = dest || new matrixArrayType(15);
        
        var eyex = eye[0],
    		eyey = eye[1],
    		eyez = eye[2],
    		upx = up[0],
    		upy = up[1],
    		upz = up[2],
    		centerx = center[0],
    		centery = center[1],
    		centerz = center[2];

    	if (eyex == centerx && eyey == centery && eyez == centerz) {
    		return mat4.identity(dest);
    	}

    	var z0,z1,z2,x0,x1,x2,y0,y1,y2,len;

    	//vec3.direction(eye, center, z);
    	z0 = eyex - center[0];
    	z1 = eyey - center[1];
    	z2 = eyez - center[2];

    	// normalize (no check needed for 0 because of early return)
    	len = 1/Math.sqrt(z0*z0 + z1*z1 + z2*z2);
    	z0 *= len;
    	z1 *= len;
    	z2 *= len;

    	//vec3.normalize(vec3.cross(up, z, x));
    	x0 = upy*z2 - upz*z1;
    	x1 = upz*z0 - upx*z2;
    	x2 = upx*z1 - upy*z0;
    	len = Math.sqrt(x0*x0 + x1*x1 + x2*x2);
    	if (!len) {
    		x0 = 0;
    		x1 = 0;
    		x2 = 0;
    	} else {
    		len = 1/len;
    		x0 *= len;
    		x1 *= len;
    		x2 *= len;
    	};

    	//vec3.normalize(vec3.cross(z, x, y));
    	y0 = z1*x2 - z2*x1;
    	y1 = z2*x0 - z0*x2;
    	y2 = z0*x1 - z1*x0;

    	len = Math.sqrt(y0*y0 + y1*y1 + y2*y2);
    	if (!len) {
    		y0 = 0;
    		y1 = 0;
    		y2 = 0;
    	} else {
    		len = 1/len;
    		y0 *= len;
    		y1 *= len;
    		y2 *= len;
    	}

    	dest[0] = x0;
    	dest[1] = y0;
    	dest[2] = z0;
    	dest[3] = 0;
    	dest[4] = x1;
    	dest[5] = y1;
    	dest[6] = z1;
    	dest[7] = 0;
    	dest[8] = x2;
    	dest[9] = y2;
    	dest[10] = z2;
    	dest[11] = 0;
    	dest[12] = -(x0*eyex + x1*eyey + x2*eyez);
    	dest[13] = -(y0*eyex + y1*eyey + y2*eyez);
    	dest[14] = -(z0*eyex + z1*eyey + z2*eyez);
    	dest[15] = 1;

    	return dest;
    },
    translate: function(mat, vec, dest){
        dest = dest || new matrixArrayType(15);
        
        var x = vec[0], y = vec[1], z = vec[2];

    	var a00 = mat[0], a01 = mat[1], a02 = mat[2], a03 = mat[3];
    	var a10 = mat[4], a11 = mat[5], a12 = mat[6], a13 = mat[7];
    	var a20 = mat[8], a21 = mat[9], a22 = mat[10], a23 = mat[11];

    	dest[0] = a00;
    	dest[1] = a01;
    	dest[2] = a02;
    	dest[3] = a03;
    	dest[4] = a10;
    	dest[5] = a11;
    	dest[6] = a12;
    	dest[7] = a13;
    	dest[8] = a20;
    	dest[9] = a21;
    	dest[10] = a22;
    	dest[11] = a23;

    	dest[12] = a00*x + a10*y + a20*z + mat[12];
    	dest[13] = a01*x + a11*y + a21*z + mat[13];
    	dest[14] = a02*x + a12*y + a22*z + mat[14];
    	dest[15] = a03*x + a13*y + a23*z + mat[15];
    	return dest;
    },
};
