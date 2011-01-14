function getShader(gl, id) {
    var shaderScript = document.getElementById(id);
    if (!shaderScript) {
        return null;
    }

    var str = "";
    var k = shaderScript.firstChild;
    while (k) {
        if (k.nodeType == 3) {
            str += k.textContent;
        }
        k = k.nextSibling;
    }

    var shader;
    if (shaderScript.type == "x-shader/x-fragment") {
        shader = gl.createShader(gl.FRAGMENT_SHADER);
    } else if (shaderScript.type == "x-shader/x-vertex") {
        shader = gl.createShader(gl.VERTEX_SHADER);
    } else {
        return null;
    }

    gl.shaderSource(shader, str);
    gl.compileShader(shader);

    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        alert(gl.getShaderInfoLog(shader));
        return null;
    }

    return shader;
}

function RenderObject() {
    this.buffer = null;
}
Object.extend(RenderObject.prototype, {
    draw: function(gl){
        gl.bindBuffer(gl.ARRAY_BUFFER, this.buffer);
        // 頂点シェーダに頂点座標データをセット
		gl.vertexAttribPointer(shaderProgram.vertexPositionAttribute, this.buffer.itemSize, gl.FLOAT, false, 0, 0);
		gl.drawArrays(gl.TRIANGLES, 0, this.buffer.numItems);
    }
});
