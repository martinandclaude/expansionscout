
// Repeat in warm, flank in cool: the same contrast REViewer uses, so the
// question "how far into the repeat does this read reach" is answered by
// colour before anything is counted.
const C={C:'#d98032',I:'#7b4b9e',P:'#c0392b',B:'#b8860b',o:'#a3a29c','.':'#d8d7d2'};
const FLANK='#6f8fae', CLIP='#b9b8b2';
const UW=8,UH=9,MH=6,GAP=4,PAD=6,MAXDRAW=110,TAIL=8,MINOMIT=10;
const BARW=660,BARH=11,FLANKPX=14,UNIT_ZOOM_MAX=80;

function rle(st,n){ // structure string -> [[char,count],...] over n units
  const out=[];
  for(let i=0;i<n;i++){
    const ch=(i<st.length)?st[i]:'o';
    if(out.length&&out[out.length-1][0]===ch) out[out.length-1][1]++;
    else out.push([ch,1]);
  }
  return out;
}

function cpgSeries(s,unitLen,n){ // -> array of {u, v} in unit coordinates
  if(!s) return [];
  const out=[];
  for(const part of s.split(',')){
    const [o,p]=part.split(':'); if(p===undefined) continue;
    const u=Number(o)/(unitLen||3);
    if(u>=0&&u<=n) out.push({u:u,v:Number(p)/100});
  }
  return out.sort((a,b)=>a.u-b.u);
}

function drawBar(rd,loc,pxPerUnit,showMeth){
  // One read as a continuous bar: flank, then the repeat run-length encoded so
  // a pure tract is one segment and an interruption is a mark inside it rather
  // than a block of its own. A clipped read gets an open end, because it says
  // "at least this long" and a closed bar would say something it does not.
  const n=Math.round(rd.n);
  const left=rd.cls!=='left_partial', right=rd.cls!=='right_partial';
  let x=PAD,parts=[];
  if(left){parts.push(`<rect x="${x}" y="0" width="${FLANKPX}" height="${BARH}" fill="${FLANK}"/>`);}
  else {parts.push(`<path d="M${x+FLANKPX} 0 L${x+2} ${BARH/2} L${x+FLANKPX} ${BARH}" fill="${CLIP}"/>`);}
  x+=FLANKPX;
  const x0=x;
  rle(rd.struct||'',n).forEach(([ch,count])=>{
    const w=Math.max(count*pxPerUnit,0.4);
    parts.push(`<rect x="${x}" y="0" width="${w}" height="${BARH}" fill="${C[ch]||C.o}"/>`);
    x+=w;
  });
  const x1=x;
  if(showMeth){
    const pts=cpgSeries(rd.cpg,loc.unit_len,n);
    pts.forEach((p,i)=>{
      const nx=x0+p.u*pxPerUnit;
      const w=Math.max((i+1<pts.length?(pts[i+1].u-p.u):1)*pxPerUnit,0.7);
      parts.push(`<rect x="${nx}" y="${BARH+1}" width="${w}" height="${MH}" fill="${methColor(p.v)}"/>`);
    });
  }
  if(right){parts.push(`<rect x="${x}" y="0" width="${FLANKPX}" height="${BARH}" fill="${FLANK}"/>`);x+=FLANKPX;}
  else {parts.push(`<path d="M${x} 0 L${x+FLANKPX-2} ${BARH/2} L${x} ${BARH}" fill="${CLIP}"/>`);x+=FLANKPX;}
  return {svg:parts.join(''),w:x,h:BARH+(showMeth?MH+1:0)};
}
const esc=s=>String(s==null?'':s).replace(/[&<>"]/g,m=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[m]));
const num=(v,d)=>v==null?'&mdash;':(d===undefined?v:Number(v).toFixed(d));
const pct=v=>v==null?'&mdash;':(100*v).toFixed(1)+'&nbsp;%';

function methColor(p){ // 0..1, blue -> neutral -> red
  const lo=[0x2a,0x78,0xd6],mid=[0xd6,0xd5,0xcf],hi=[0xd8,0x3a,0x35];
  let c; if(p<0.5){const t=p/0.5;c=lo.map((v,i)=>Math.round(v+(mid[i]-v)*t));}
  else{const t=(p-0.5)/0.5;c=mid.map((v,i)=>Math.round(v+(hi[i]-v)*t));}
  return '#'+c.map(v=>v.toString(16).padStart(2,'0')).join('');
}

function parseCpg(s,unitLen){ // "off:pct,..." -> {unitIndex:[pcts]}
  const out={}; if(!s) return out;
  for(const part of s.split(',')){
    const [o,p]=part.split(':'); if(p===undefined) continue;
    const u=Math.floor(Number(o)/(unitLen||3));
    (out[u]=out[u]||[]).push(Number(p));
  }
  const m={}; for(const k in out) m[k]=out[k].reduce((a,b)=>a+b,0)/out[k].length/100;
  return m;
}

function drawRead(rd,loc,showMeth,cap){
  const n=Math.round(rd.n), st=rd.struct||'';
  const partial=rd.cls==='left_partial'||rd.cls==='right_partial';
  let head=n,tail=0,omit=0,brk=-1;
  if(n>cap&&(n-cap)>=MINOMIT){head=cap-TAIL;tail=TAIL;omit=n-head-tail;brk=head;}
  const idx=[]; for(let i=0;i<head;i++)idx.push(i);
  for(let i=n-tail;i<n;i++)idx.push(i);
  const cpg=showMeth?parseCpg(rd.cpg,loc.unit_len):null;
  let x=PAD,parts=[],tips=[];
  idx.forEach((i,k)=>{
    if(brk>=0&&k===brk){
      parts.push(`<path d="M${x+2} -1 L${x+7} ${UH+1} M${x+7} -1 L${x+12} ${UH+1}" stroke="#a3a29c" stroke-width="1.2" fill="none"/>`);
      parts.push(`<text x="${x+15}" y="${UH-1}" font-size="9.5" fill="#6b6a66">${omit} more</text>`);
      x+=20+String(omit).length*6+34;
    }
    const ch=(i<st.length)?st[i]:'o';
    parts.push(`<rect x="${x}" y="0" width="${UW-1}" height="${UH}" rx="1.5" fill="${C[ch]||C.o}"/>`);
    if(cpg&&cpg[i]!==undefined)
      parts.push(`<rect x="${x}" y="${UH+1}" width="${UW-1}" height="${MH}" rx="1" fill="${methColor(cpg[i])}"/>`);
    x+=UW;
  });
  if(partial){ // an open-ended bar: this read says "at least this long"
    parts.push(`<path d="M${x+1} 0 L${x+7} ${UH/2} L${x+1} ${UH}" fill="none" stroke="#6b6a66" stroke-width="1.2"/>`);
    x+=12;
  }
  const h=UH+(showMeth?MH+1:0);
  const label=(rd.hp?`H${rd.hp}`:'')+(partial?' ≥':'');
  return {svg:`<g transform="translate(0,0)">${parts.join('')}</g>`,w:x,h:h,label:label};
}

function pileup(loc,opts){
  const showMeth=opts.meth&&loc.reads.some(r=>r.cpg);
  let reads=loc.reads.slice();
  if(opts.sort==='meth') reads.sort((a,b)=>(b.meth==null?-1:b.meth)-(a.meth==null?-1:a.meth));
  else reads.sort((a,b)=>b.n-a.n||(a.id<b.id?-1:1));
  const maxN=Math.max(1,...reads.map(r=>r.n));
  // One scale for the whole locus, so allele lengths are comparable by eye.
  // A 17 next to a 700 is a sliver, which is the honest picture and the one a
  // Southern blot gives.
  const pxPerUnit=Math.min(BARW/maxN,6);
  const groups={};
  reads.forEach(r=>{(groups[r.allele||0]=groups[r.allele||0]||[]).push(r)});
  // Unassigned reads last: they are the leftovers, not the headline. Putting
  // them first buried allele 1 below a group nobody asked about.
  const keys=Object.keys(groups).sort((a,b)=>(Number(a)||99)-(Number(b)||99));
  // One canvas width for every group, so the unit counts line up in a column
  // down the locus instead of chasing each group's longest bar.
  const canvas=Math.max(...reads.map(r=>(opts.view==='units'
      ? drawRead(r,loc,false,opts.full?1e6:MAXDRAW).w
      : drawBar(r,loc,pxPerUnit,false).w)))+74;
  let out=[];
  keys.forEach(k=>{
    const al=(loc.alleles||[]).find(a=>a.n==='a'+k);
    let head='';
    if(al){
      const cls=al.cls||'';
      head=`<div style="margin:14px 0 5px"><b>allele ${k}</b> &nbsp;`+
        `<span class="cls ${esc(cls)}">${esc(cls.replace('_',' '))}</span> &nbsp;`+
        `<span class="muted">${num(al.median,0)} units, ${num(al.support,0)} reads`+
        (al.p5!=null?`, reads span ${num(al.p5,0)}–${num(al.p95,0)}`:'')+
        (al.meth_up!=null?`, promoter 5mC ${pct(al.meth_up)}`:
          (al.meth_tract!=null?`, tract 5mC ${pct(al.meth_tract)}`:''))+
        `</span></div>`;
      if(al.dispersed) head+=`<div class="why" style="margin:0 0 6px">These reads were not `+
        `summarised as one allele: the spread is too wide for a single length to describe `+
        `them, so no clinical band was assigned. The range above is the result.</div>`;
    } else {
      const np=groups[k].filter(r=>r.n_path>0);
      head=`<div style="margin:14px 0 5px"><b>unassigned reads</b> `+
        `<span class="muted">${groups[k].length} reads that did not join either `+
        `allele: too few to call, or clipped</span></div>`;
      if(np.length) head+=`<div class="why" style="margin:0 0 6px">`+
        `${np.length===1?'One of them carries':np.length+' of them carry'} `+
        `pathogenic-motif units `+
        `(${np.map(r=>Math.round(r.n_path)).join(', ')}). Below the support `+
        `threshold, so they are in no allele and in no count -- shown because `+
        `a table cannot show them.</div>`;
    }
    const rows=groups[k].map(r=>({
      d: opts.view==='units' ? drawRead(r,loc,showMeth,opts.full?1e6:MAXDRAW)
                             : drawBar(r,loc,pxPerUnit,showMeth),
      r: r}));
    const w=canvas;
    let y=0,body='';
    rows.forEach(x=>{
      body+=`<g transform="translate(56,${y})" class="rd" data-id="${esc(x.r.id)}" `+
        `data-n="${x.r.n}" data-cls="${esc(x.r.cls)}" data-hp="${x.r.hp}" `+
        `data-mapq="${x.r.mapq==null?'':x.r.mapq}" data-meth="${x.r.meth==null?'':x.r.meth}" `+
        `data-lb="${x.r.lb==null?'':x.r.lb}" `+
        `data-npath="${x.r.n_path==null?'':x.r.n_path}">${x.d.svg}</g>`;
      const lbl=(x.r.hp?'H'+x.r.hp:'');
      body+=`<text x="0" y="${y+BARH-2}" font-size="10" fill="#6b6a66">${esc(lbl)}</text>`;
      const partial=x.r.cls==='left_partial'||x.r.cls==='right_partial';
      body+=`<text x="${w-12}" y="${y+BARH-2}" font-size="10" fill="#6b6a66" `+
        `text-anchor="end">${partial?'&#8805;':''}${Math.round(x.r.n)}</text>`;
      y+=x.d.h+GAP;
    });
    out.push(head+`<div class="pileup"><svg width="${w}" height="${y}" `+
      `xmlns="http://www.w3.org/2000/svg">${body}</svg></div>`);
  });
  return out.join('');
}

function bandTable(loc){
  const bs=loc.bands||[]; if(!bs.length) return '<span class="muted">no bands published for this locus</span>';
  const cite=k=>{const s=(DATA.sources||{})[k]; if(!s) return esc(k||'');
    return `<span title="${esc(s.citation||'')}">${esc(s.citation||k)}</span>`+
      (s.jurisdiction?` <span class="muted">(${esc(s.jurisdiction)})</span>`:'')+
      (s.pmcid?` <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/${esc(s.pmcid)}/">${esc(s.pmcid)}</a>`:'');};
  let prev=null,rows='';
  bs.forEach(b=>{
    // A band that states its own lower bound is a published range; one that
    // gives only an upper bound is part of an ordered partition and starts
    // where the previous band ended.
    const lo=(b.lower!=null)?b.lower:(prev==null?0:prev+1), hi=b.upper;
    if(b.lower!=null&&prev!=null&&b.lower>prev+1)
      rows+=`<tr class="gap"><td class="muted">&mdash;</td><td class="muted">`+
        `${prev+1}${b.lower-1>prev+1?'–'+(b.lower-1):''}</td>`+
        `<td class="muted">no band published for this range</td></tr>`;
    const rng=hi==null?`${lo} and above`:(lo===hi?`${lo}`:`${lo}–${hi}`);
    const q=b.boundary_quote||b.naming_quote;
    rows+=`<tr><td>${esc(b.label||'')}${b.inhouse?' <span class="tag">in-house</span>':''}</td>`+
      `<td>${rng}</td><td class="muted">${cite(b.boundary)}`+
      (b.naming&&b.naming!==b.boundary?`<br>named by ${cite(b.naming)}`:'')+
      (q?`<br><i>&ldquo;${esc(q)}&rdquo;</i>`:'')+`</td></tr>`;
    prev=(hi==null)?prev:hi;
  });
  return `<table class="bands"><tr><th>band</th><th>units</th><th>where it comes from</th></tr>${rows}</table>`;
}

function detailPanel(loc){
  const ev={NONE:'no partial-read evidence',PARTIALS_PRESENT:'partial reads present but not conclusive',
            EXPANSION_LB:'expansion supported by partial reads'};
  const kv=[
    ['reads',`${num(loc.n_spanning,0)} spanning, ${num(loc.n_split,0)} split, ${num(loc.n_partial,0)} partial`],
    ['motif',`${esc(loc.motif||'')} on the gene strand, ${esc(loc.ref_motif||'')} on the reference`],
    ['assignment',esc(loc.method||'')+(loc.tagged_fraction!=null?` &middot; ${pct(loc.tagged_fraction)} of reads carry a haplotag`:'')],
  ];
  if(loc.lb!=null) kv.push(['lower bound',`at least ${num(loc.lb,0)} units, from ${num(loc.lb_support,0)} reads`]);
  if(loc.evidence&&loc.evidence!=='NONE') kv.push(['evidence',esc(ev[loc.evidence]||loc.evidence)]);
  if(loc.detect!=null) kv.push(['detectability',
    `${pct(loc.detect)} of reads are long enough to span a pathogenic allele here`]);
  if(loc.negative_reliable===0) kv.push(['caution',
    'a normal call at this locus is not evidence of a normal result']);
  if(loc.notes) kv.push(['notes',esc(loc.notes)]);
  const maxN=Math.max(1,...loc.reads.map(r=>r.n));
  const zoomable=maxN<=UNIT_ZOOM_MAX;
  const ctl=`<div class="controls">
     <span><label>view </label><select data-ctl="view"><option value="bars">continuous</option>
       <option value="units"${zoomable?'':' disabled'}>unit by unit</option></select>
       ${zoomable?'':`<span class="muted"> &mdash; unit view is off above `+
         `${UNIT_ZOOM_MAX} units, where it stops being readable</span>`}</span>
     <span><label>sort </label><select data-ctl="sort"><option value="length">by length</option>
       <option value="meth">by methylation</option></select></span>
     <span><label>methylation </label><select data-ctl="meth"><option value="on">show</option>
       <option value="off">hide</option></select></span></div>`;
  const leg=`<div class="legend">
     <span><i class="sw" style="background:${FLANK}"></i>flank (schematic, not to scale)</span>
     <span><i class="sw" style="background:${C.C}"></i>repeat</span>
     <span><i class="sw" style="background:${C.I}"></i>interruption</span>
     <span><i class="sw" style="background:${C.P}"></i>pathogenic motif</span>
     <span><i class="sw" style="background:${C.o}"></i>other</span>
     <span><i class="sw" style="background:${methColor(0)}"></i>unmethylated</span>
     <span><i class="sw" style="background:${methColor(1)}"></i>methylated</span>
     <span><i class="sw" style="background:${CLIP}"></i>clipped end: length is a lower bound</span></div>`;
  return `<div class="panel"><h3>${esc(loc.gene)} &middot; ${esc(loc.id)}</h3>
    <p class="why">One read per row, longest first, on a single scale for this
    locus so allele lengths compare by eye. Flanks are drawn as fixed caps to
    show whether a read spans the repeat; they are not to scale, because on
    long reads they would be kilobases and would swamp the tract.</p>
    ${ctl}<div data-pileup></div>${leg}
    <dl class="kv">${kv.map(([k,v])=>`<dt>${k}</dt><dd>${v}</dd>`).join('')}</dl>
    <h3 style="margin-top:16px">Clinical bands, and where each one comes from</h3>
    <p class="why">The boundary and the name of a band can come from different
    places -- a number from the catalogue with a name from a guideline is the
    common case -- so both are given. Guidelines differ between countries.</p>
    ${bandTable(loc)}</div>`;
}

document.addEventListener('DOMContentLoaded',function(){
  const byId={}; DATA.loci.forEach(l=>byId[l.id]=l);
  document.querySelectorAll('tr.locus').forEach(tr=>{
    tr.addEventListener('click',function(){
      const id=tr.dataset.locus, det=tr.nextElementSibling;
      if(det.dataset.built!=='1'){
        const loc=byId[id];
        det.querySelector('td').innerHTML=detailPanel(loc);
        const host=det.querySelector('[data-pileup]');
        const opts={sort:'length',meth:true,full:false,view:'bars'};
        const redraw=()=>{host.innerHTML=pileup(loc,opts)};
        det.querySelectorAll('[data-ctl]').forEach(sel=>{
          sel.addEventListener('change',e=>{
            const k=sel.dataset.ctl;
            opts[k]=(k==='sort'||k==='view')?sel.value:(sel.value==='on');
            if(k==='meth') opts.meth=(sel.value==='on');
            redraw(); e.stopPropagation();
          });
          sel.addEventListener('click',e=>e.stopPropagation());
        });
        redraw();
        det.dataset.built='1';
      }
      det.hidden=!det.hidden;
    });
  });
  document.querySelectorAll('.coord').forEach(b=>{
    b.addEventListener('click',function(e){
      e.stopPropagation();
      navigator.clipboard&&navigator.clipboard.writeText(b.dataset.coord);
      const t=b.textContent; b.textContent='copied'; setTimeout(()=>b.textContent=t,800);
    });
  });
  const tip=document.getElementById('tip');
  document.addEventListener('mouseover',function(e){
    const g=e.target.closest&&e.target.closest('.rd'); if(!g){tip.style.display='none';return;}
    const d=g.dataset;
    tip.innerHTML=`<b>${esc(d.id)}</b><br>${Math.round(d.n)} units &middot; ${esc(d.cls)}`+
      (d.hp&&d.hp!=='0'?` &middot; haplotype ${esc(d.hp)}`:'')+
      (d.mapq?`<br>MAPQ ${esc(d.mapq)}`:'')+
      (d.meth?`<br>tract 5mC ${(100*Number(d.meth)).toFixed(1)} %`:'')+
      (d.lb?`<br>lower bound ${esc(d.lb)} units`:'')+
      (d.npath&&Number(d.npath)>0?`<br><b>${esc(d.npath)} pathogenic-motif units</b>`:'');
    tip.style.display='block';
  });
  document.addEventListener('mousemove',function(e){
    tip.style.left=Math.min(e.clientX+14,window.innerWidth-360)+'px';
    tip.style.top=(e.clientY+16)+'px';
  });
});
