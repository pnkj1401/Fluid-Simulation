export default class Smoke{

    constructor(dt, diffusion, viscosity){
        this.size = N;
        this.dt = dt;
        this.diff = diffusion;
        this.visc = viscosity;

        this.s = new Array(N*N).fill(0);
        this.density = new Array(N*N).fill(0);
        
        this.Vx = new Array(N*N).fill(0);
        this.Vy = new Array(N*N).fill(0);

        this.Vx0 = new Array(N*N).fill(0);
        this.Vy0 = new Array(N*N).fill(0);
 
    }

    step(){
        let visc = this.visc;
        let diff = this.diff;
        let dt = this.dt;
        let Vx = this.Vx;
        let Vy = this.Vy;
        let Vx0 = this.Vx0;
        let Vy0 = this.Vy0;
        let s = this.s;
        let density = this.density;
        
        
        diffuse(1, Vx0, Vx, visc, dt);
        diffuse(2, Vy0, Vy, visc, dt);
        
        project(Vx0, Vy0, Vx, Vy);
        
        advect(1, Vx, Vx0, Vx0, Vy0, dt);
        advect(2, Vy, Vy0, Vx0, Vy0, dt);
        
        project(Vx, Vy, Vx0, Vy0);
        
        diffuse(0, s, density, diff, dt);
        advect(0, density, s, Vx, Vy, dt);
        
    }

    addDensity(x, y, amount){
      const index = IX(x, y);
      if(index > N*N && this.density[index] > 1024 ) return;
      this.density[index] = this.density[index] + amount;
      
    }
    
    addVelocity(x, y, amountX, amountY){
        const index = IX(x,y);
        this.Vx[index] += amountX;
        this.Vy[index] += amountY;
    }
    
    render(){ 
        for(let i=0; i<N; i++){
            for(let j=0; j<N; j++){
                let x = i * SCALE + i * 0;
                let y = j * SCALE + j * 0;
                let d = this.density[IX(i, j)];
                const alpha = Math.floor(d / 1024 * 256).toString(16).padStart(2, "0");
                ctx.fillStyle = "#0dd0ff" + alpha;
                ctx.fillRect(x, y, SCALE, SCALE);
            }
        }
    }

}

function IX(x, y){
    x = Math.max(0,Math.min(x, N-1));
    y = Math.max(0,Math.min(y, N-1));
    return Math.floor(x + y * N);
}


function diffuse(b, x, x0, diff, dt){
     let a = dt * diff * (N-2) * (N-2);
     linSolve(b, x, x0, a, 1 + 6 *a);
}

function linSolve(b, x, x0, a, c){
     let cRecip = 1.0 / c;
     for(let k=0; k<iter; k++){
        for(let j=1; j<N-1; j++){
            for(let i=1; i<N-1; i++){
                const index = IX(i, j);
                 x[index] = (
                    x0[IX(i, j)] + 
                    a * (x[IX(i+1, j)]
                        + x[IX(i-1, j)] 
                        + x[IX(i, j+1)]
                        + x[IX(i, j-1)]
                        )
                 ) * cRecip;
            }
        }
        set_bnd(b, x);
     }
}

function project(velocX, velocY, p, div){
    for(let j=1;j<N-1;j++){
        for(let i=1;i<N-1;i++){
            div[IX(i,j)] = -0.5 * (
                velocX[IX(i+1, j)]
                -velocX[IX(i-1, j)]
                +velocY[IX(i, j+1)]
                -velocY[IX(i, j-1)]) / N;
            p[IX(i,j)] = 0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);
    linSolve(0, p, div, 1, 6);

    for(let j=1; j<N-1; j++){
        for(let i=1; i<N-1; i++){
            velocX[IX(i, j)] -= 0.5 * (p[IX(i+1, j)] - p[IX(i-1, j)]) * N;
            velocY[IX(i, j)] -= 0.5 * (p[IX(i, j+1)] - p[IX(i, j-1)]) * N
        }
    }
    set_bnd(1, velocX);
    set_bnd(2, velocY);
}

function advect(b, d, d0, velocX, velocY, dt)
{
    let i0, i1, j0, j1;
    
    let dtx = dt * (N - 2);
    let dty = dt * (N - 2);
    
    let s0, s1, t0, t1;
    let tmp1, tmp2, x, y;
    
    let Nfloat = N;
    let ifloat, jfloat;
    let i, j;
    
        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[IX(i, j)];
                tmp2 = dty * velocY[IX(i, j)];
                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                
                if(x < 0.5) x = 0.5; 
                if(x > Nfloat + 0.5) x = Nfloat + 0.5; 
                i0 = Math.floor(x); 
                i1 = i0 + 1.0;
                if(y < 0.5) y = 0.5; 
                if(y > Nfloat + 0.5) y = Nfloat + 0.5; 
                j0 = Math.floor(y);
                j1 = j0 + 1.0; 
               
                
                s1 = x - i0; 
                s0 = 1.0 - s1; 
                t1 = y - j0; 
                t0 = 1.0 - t1;
                
                let i0i = Number.parseInt(i0);
                let i1i = Number.parseInt(i1);
                let j0i = Number.parseInt(j0);
                let j1i = Number.parseInt(j1);
                const _0 = d0[IX(i0i, j0i)];
                const _1 = (d0[IX(i0i, j1i)]);
                const _2 = (d0[IX(i1i, j0i)]);
                const _3 = (d0[IX(i1i, j1i)]);
                let val = s0 * (
                    t0 * _0
                + t1 * _1)
                 +s1 * (
                    t0 * _2
                +( t1 * _3)
                );
                if(Number.isNaN(val)) debugger;
                d[IX(i, j)] = val;
                
                    
            }
        }

    set_bnd(b, d);
}

function set_bnd(b, x)
{
    for(let i = 1; i < N - 1; i++) {
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N-1)] = b == 2 ? -x[IX(i, N-2)] : x[IX(i, N-2)];
    }

    for(let j = 1; j < N - 1; j++) {
        x[IX(0  , j)] = b == 1 ? -x[IX(1  , j)] : x[IX(1  , j)];
        x[IX(N-1, j)] = b == 1 ? -x[IX(N-2, j)] : x[IX(N-2, j)];
    }
    
    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N-1)] = 0.5 * (x[IX(1, N-1)] + x[IX(0, N-2)]);
    x[IX(N-1, 0)] = 0.5 * (x[IX(N-2, 0)] + x[IX(N-1, 1)]);
    x[IX(N-1, N-1)] = 0.5 * (x[IX(N-2, N-1)] + x[IX(N-1, N-2)]);

}