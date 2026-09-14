/* influence-lines.js — 影響線／設計效應共用核心（index.html 影響線計算器 ＋ analyzer.html ⑦ 連續梁）
 * 數值剛度法：支援簡支與連續梁；HS20-44 軸組雙向掃描；剪力影響線於斷面存左右極限（跳躍 1）。
 * 簡支結果對 Python bridgecalc（influence.py）交叉驗證——見 index.html 的自驗證面板。
 * 掛 window.IL；node 則 module.exports。
 */
(function (global) {
  'use strict';
 var TRK={P:[36,144,144],x:[0,4.25,8.5]}, LANE=9.4, PM=80, PV=116;
 // 車可雙向行駛：TRKS[0]＝原向、[1]＝掉頭（軸序與軸距鏡射）；單向會使非跨中彎矩與短跨剪力偏低
 var TRKS=[TRK,{P:TRK.P.slice().reverse(),x:TRK.x.map(function(d){return TRK.x[2]-d}).reverse()}];
 // ── 簡支解析影響線（與 engine/influence.py 一致，用於自驗證）──
 function ilPeak(L,a){ return a*(L-a)/L; }
 function ilMa(L,a,p){ return p<=a ? p*(L-a)/L : a*(L-p)/L; }
 // 台灣 HS20-44 設計卡車（自驗證用，與 engine/influence.py 的 taiwan_* 一致）
 var TWa=[36,144,144], TWs=[0,4.25,8.5];
 function absMaxTW(L){ var best=0; for(var a=0.2;a<L;a+=0.2){ for(var s=-TWs[2];s<=L;s+=0.2){ var t=0; for(var i=0;i<3;i++){ var x=s+TWs[i]; if(x>=0&&x<=L)t+=TWa[i]*ilMa(L,a,x);} if(Math.abs(t)>Math.abs(best))best=t; } } return best; }
 function twLaneM(L){ return 9.4*L*L/8 + 80*L/4; }              // 車道均布 + 集中載重(彎矩)
 function twImpact(L){ return Math.min(0.30, 15.24/(L+38.1)); }
 function twPerLane(L){ return Math.max(absMaxTW(L), twLaneM(L))*(1+twImpact(L)); }   // 取大×(1+I)
 // ── 數值剛度法影響線（連續梁；簡支亦適用）──
 function inv(A){var n=A.length,M=A.map(function(r,i){return r.concat(Array.from({length:n},function(_,j){return i===j?1:0}))});for(var i=0;i<n;i++){var p=i;for(var r=i+1;r<n;r++)if(Math.abs(M[r][i])>Math.abs(M[p][i]))p=r;var t=M[i];M[i]=M[p];M[p]=t;var pv=M[i][i];for(var j=0;j<2*n;j++)M[i][j]/=pv;for(var r2=0;r2<n;r2++)if(r2!==i){var f=M[r2][i];for(var j2=0;j2<2*n;j2++)M[r2][j2]-=f*M[i][j2]}}return M.map(function(r){return r.slice(n)})}
 // 每跨分割段數：包絡只在節點取值，密度不足會低估峰值（簡支 24 段時較絕對最大低 0.1%）
 function model(spans,neIn){var ne=neIn||(spans.length>1?24:60),nodes=[0],sx=[0],cum=0;spans.forEach(function(Ln){for(var k=1;k<=ne;k++)nodes.push(+(cum+Ln*k/ne).toFixed(4));cum+=Ln;sx.push(+cum.toFixed(4))});var sIdx=sx.map(function(x){var b=0;for(var i=0;i<nodes.length;i++)if(Math.abs(nodes[i]-x)<Math.abs(nodes[b]-x))b=i;return b});return {nodes:nodes,sx:sx,sIdx:sIdx,tot:cum}}
 // rs=true：支承節點取「右側」剪力（計入該支承反力）；預設為左側。僅對 eff=V 有意義
 // 剛度矩陣與其反矩陣只與模型有關 → 每個模型算一次並快取（原本每個斷面都重算一次反矩陣）
 function prepare(m){
  if(m._ctx) return m._ctx;
  var nodes=m.nodes,N=nodes.length,nd=2*N,K=[];
  for(var i=0;i<nd;i++)K.push(new Array(nd).fill(0));
  for(var e=0;e<N-1;e++){
   var le=nodes[e+1]-nodes[e],cc=1/(le*le*le),
       ke=[[12*cc,6*le*cc,-12*cc,6*le*cc],[6*le*cc,4*le*le*cc,-6*le*cc,2*le*le*cc],
           [-12*cc,-6*le*cc,12*cc,-6*le*cc],[6*le*cc,2*le*le*cc,-6*le*cc,4*le*le*cc]],
       mp=[2*e,2*e+1,2*e+2,2*e+3];
   for(var a=0;a<4;a++)for(var b=0;b<4;b++)K[mp[a]][mp[b]]+=ke[a][b];
  }
  var cons={}; m.sIdx.forEach(function(s){cons[2*s]=1});
  var free=[],idxOf=new Array(nd).fill(-1);
  for(var d=0;d<nd;d++)if(!cons[d]){idxOf[d]=free.length;free.push(d)}
  var Kr=free.map(function(a){return free.map(function(b){return K[a][b]})});
  m._ctx={N:N,nd:nd,K:K,cons:cons,free:free,idxOf:idxOf,Kri:inv(Kr)};
  return m._ctx;
 }
 // rs=true：支承節點取「右側」剪力（計入該支承反力）；預設為左側。僅對 eff=V 有意義
 function solveIL(m,tIdx,eff,rs){
  var c=prepare(m),nodes=m.nodes,N=c.N,nd=c.nd,K=c.K,free=c.free,Kri=c.Kri;
  var xt=nodes[tIdx],il=[];
  for(var p=0;p<N;p++){
   if(c.cons[2*p]){il.push(0);continue}
   // 單位載重的解 = −Kri 的對應行（原本做一次矩陣×單位向量，等價但慢）
   var col=c.idxOf[2*p],dfull=new Array(nd).fill(0);
   for(var k=0;k<free.length;k++)dfull[free[k]]=-Kri[k][col];
   var R=m.sIdx.map(function(s){var sm=0;for(var j=0;j<nd;j++)sm+=K[2*s][j]*dfull[j];return sm});
   var val=0,px=nodes[p];
   if(eff==='M'){
    for(var k2=0;k2<m.sx.length;k2++)if(m.sx[k2]<xt-1e-6)val+=R[k2]*(xt-m.sx[k2]);
    if(px<xt-1e-6)val-=(xt-px);
   } else {
    var lim=rs?xt+1e-6:xt-1e-6;
    for(var k3=0;k3<m.sx.length;k3++)if(m.sx[k3]<lim)val+=R[k3];
    if(px<lim)val-=1;
   }
   il.push(val);
  }
  var ref=(eff==='M'?Math.max(m.tot,1)/4:1),out=il.map(function(v){return Math.abs(v)<ref*1e-3?0:v});
  // 剪力影響線於斷面有跳躍 1：以重複節點存左/右極限（xs/vs），否則內插會把跳躍抹進鄰元素、支承旁的軸被低估
  if(eff==='V'){
   var sup=m.sIdx.indexOf(tIdx)>=0,jl,jr;
   if(sup){jl=rs?0:-1;jr=rs?1:0}else{jr=out[tIdx];jl=jr-1}
   if(tIdx===0)jl=jr; if(tIdx===N-1)jr=jl;   // 梁端外側不存在
   out.xs=nodes.slice(0,tIdx+1).concat(nodes.slice(tIdx));
   out.vs=out.slice(0,tIdx).concat([jl,jr],out.slice(tIdx+1));
  }
  return out;
 }
 function interp(nodes,il,x){if(x<=nodes[0])return il[0];if(x>=nodes[nodes.length-1])return il[il.length-1];for(var i=0;i<nodes.length-1;i++)if(x>=nodes[i]&&x<=nodes[i+1]){var t=(x-nodes[i])/(nodes[i+1]-nodes[i]);return il[i]*(1-t)+il[i+1]*t}return 0}
 // 卡車軸組於 s（最左軸位置）、車向 dir 時之效應；軸位夾到 [0,tot] 免支承上的軸被浮點誤差略去
 // 回傳陣列：剪力影響線有跳躍時於 s±δ 各算一次（軸恰在斷面時左右極限都取到）
 function truckSum(nodes,il,tot,s,dir){var T=TRKS[dir||0],X=il.xs||nodes,Y=il.vs||il,ds=il.xs?[-1e-7,1e-7]:[0];return ds.map(function(e){var sum=0;for(var a=0;a<3;a++){var ax=s+e+T.x[a];if(ax>=-1e-6&&ax<=tot+1e-6)sum+=T.P[a]*interp(X,Y,Math.min(Math.max(ax,0),tot))}return sum})}
 function truckScan(tot,step,fn){var span=TRK.x[2],n=Math.round((tot+span)/step);for(var d=0;d<2;d++)for(var k=0;k<=n;k++)fn(-span+k*step,d)}
 function absMaxOf(a){return a.reduce(function(b,v){return Math.abs(v)>Math.abs(b)?v:b},0)}
 function truckMax(nodes,il,tot){var best=0,bx=0,bd=0;truckScan(tot,0.25,function(s,d){var sum=absMaxOf(truckSum(nodes,il,tot,s,d));if(Math.abs(sum)>Math.abs(best)){best=sum;bx=s;bd=d}});return {v:best,x:bx,dir:bd}}
 function laneEnv(nodes,il,P){var X=il.xs||nodes,Y=il.vs||il,posA=0,negA=0,mx=0,mn=0;for(var i=0;i<X.length-1;i++){var dx=X[i+1]-X[i],av=(Y[i]+Y[i+1])/2;if(av>0)posA+=av*dx;else negA+=av*dx}Y.forEach(function(v){if(v>mx)mx=v;if(v<mn)mn=v});return {pos:LANE*posA+P*mx,neg:LANE*negA+P*mn}}
 function laneMax(nodes,il,P){var e=laneEnv(nodes,il,P);return Math.abs(e.pos)>=Math.abs(e.neg)?e.pos:e.neg}
 var lf={'1':[1,1.0],'2':[2,1.0],'2b':[2,1.0],'4':[4,0.75]};
 // ── 包絡線：每斷面的最大正/負設計效應 ──
 function truckEnv(nodes,il,tot){var pos=0,neg=0;truckScan(tot,0.25,function(s,d){truckSum(nodes,il,tot,s,d).forEach(function(sum){if(sum>pos)pos=sum;if(sum<neg)neg=sum})});return {pos:pos,neg:neg}}
 function isSup(m,t){return m.sIdx.indexOf(t)>=0}
 // 支承節點剪力有左右兩側（右側計入該支承反力），兩側都算
 function sideILs(m,t,eff){var a=[solveIL(m,t,eff)];if(eff==='V'&&isSup(m,t))a.push(solveIL(m,t,eff,true));return a}
 function envelope(m,eff,I,nl,red){var P=eff==='M'?PM:PV,fac=(1+I)*nl*red,posE=[],negE=[];for(var t=0;t<m.nodes.length;t++){var pos=0,neg=0;sideILs(m,t,eff).forEach(function(il){var te=truckEnv(m.nodes,il,m.tot),le=laneEnv(m.nodes,il,P);pos=Math.max(pos,te.pos,le.pos);neg=Math.min(neg,te.neg,le.neg)});posE.push(pos*fac);negE.push(neg*fac)}return {pos:posE,neg:negE}}

  // 均布載重效應 = w × ∫η dx（梯形積分；同時給正/負面積供包絡用）
  function ilArea(nodes, il) {
    var X = il.xs || nodes, Y = il.vs || il, pos = 0, neg = 0;
    for (var i = 0; i < X.length - 1; i++) {
      var dx = X[i + 1] - X[i], av = (Y[i] + Y[i + 1]) / 2;
      if (av > 0) pos += av * dx; else neg += av * dx;
    }
    return { total: pos + neg, pos: pos, neg: neg };
  }
  // 區段均布載重效應 = w × ∫[x1,x2] η dx（梯形逐段裁切；等效載重法算 M_total 用）
  function ilAreaRange(nodes, il, x1, x2) {
    var X = il.xs || nodes, Y = il.vs || il, s = 0;
    for (var i = 0; i < X.length - 1; i++) {
      var a = X[i], b = X[i + 1]; if (b <= a) continue;
      var lo = Math.max(a, x1), hi = Math.min(b, x2); if (hi <= lo) continue;
      var t1 = (lo - a) / (b - a), t2 = (hi - a) / (b - a);
      s += ((Y[i] + (Y[i + 1] - Y[i]) * t1) + (Y[i] + (Y[i + 1] - Y[i]) * t2)) / 2 * (hi - lo);
    }
    return s;
  }
  var IL = { prepare: prepare, ilAreaRange: ilAreaRange, TRK: TRK, TRKS: TRKS, LANE: LANE, PM: PM, PV: PV,
             inv: inv, model: model, solveIL: solveIL, interp: interp,
             truckSum: truckSum, truckScan: truckScan, absMaxOf: absMaxOf,
             truckMax: truckMax, truckEnv: truckEnv, laneEnv: laneEnv, laneMax: laneMax,
             isSup: isSup, sideILs: sideILs, envelope: envelope, ilArea: ilArea,
             ilPeak: ilPeak, ilMa: ilMa, absMaxTW: absMaxTW, twLaneM: twLaneM,
             twImpact: twImpact, twPerLane: twPerLane };
  global.IL = IL;
  if (typeof module !== 'undefined' && module.exports) module.exports = IL;
})(typeof self !== 'undefined' ? self : this);
