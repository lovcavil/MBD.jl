using Cassette: @context ,@overdub
Cassette.@context PrintCtx
Cassette.prehook(::PrintCtx, f, args...) = println(f, args)
Cassette.overdub(PrintCtx(), /, 1, 2)