import Smoke from "./smoke.js";
/** @type {HTMLCanvasElement} */
const canvas = document.getElementById('canvas');

globalThis.ctx = canvas.getContext("2d");
globalThis.N = 150;
globalThis.iter = 4;


let isClicked = false;

canvas.width = 800;
canvas.height = 800;
globalThis.SCALE = canvas.width / N;
let fluid = new Smoke(0.05, 0, 0);    

let mouse = { x:0, y: 0 };
canvas.addEventListener("mousedown",(event)=>{
    isClicked = true;
    mouse.x = Math.floor(event.offsetX/SCALE);
    mouse.y = Math.floor(event.offsetY/SCALE);
})
canvas.addEventListener("mouseup",(event)=>{
    isClicked = false;
})
canvas.addEventListener("mousemove",(event)=>{
    mouse.x = Math.floor(event.offsetX/SCALE);
    mouse.y = Math.floor(event.offsetY/SCALE);
})
draw();   

function draw(){
    ctx.clearRect(0, 0, canvas.width, canvas.height)
    fluid.step();
    fluid.render();
    if(isClicked){
        fluid.addDensity(mouse.x, mouse.y, 1000);
        fluid.addVelocity(mouse.x, mouse.y, 0.25, 0.25);
    }
    requestAnimationFrame(draw);
}