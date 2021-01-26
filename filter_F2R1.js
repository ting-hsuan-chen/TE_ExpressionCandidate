{
"filters" : [    
  {"id" : "FirstinPair",
   "isFirstMate" : "true"
  },     
  {"id" : "F",
   "isReverseStrand" : "false"
  },
  {"id" : "R",
   "isReverseStrand" : "true"
  },
  {"id" : "SecondinPair",
   "isFirstMate" : "false"
  },
  {"id" : "P",
   "isPaired" : "true"
  }
],
"rule": "( P & FirstinPair & R) | ( P & SecondinPair & F)"
}