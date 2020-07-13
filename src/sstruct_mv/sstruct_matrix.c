/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision$
 ***********************************************************************EHEADER*/

/******************************************************************************
 *
 * Member functions for hypre_SStructPMatrix class.
 *
 *****************************************************************************/

#include "_hypre_sstruct_mv.h"

/*==========================================================================
 * SStructPMatrix routines
 *==========================================================================*/

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructPMatrixRef( hypre_SStructPMatrix  *matrix,
                         hypre_SStructPMatrix **matrix_ref )
{
   hypre_SStructPMatrixRefCount(matrix) ++;
   *matrix_ref = matrix;

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructPMatrixCreate( MPI_Comm               comm,
                            hypre_SStructPGrid    *pgrid,
                            hypre_SStructStencil **stencils,
                            hypre_SStructPMatrix **pmatrix_ptr )
{
   HYPRE_Int              ndim  = hypre_SStructPGridNDim(pgrid);
   HYPRE_Int              nvars = hypre_SStructPGridNVars(pgrid);

   hypre_SStructPMatrix  *pmatrix;
   HYPRE_Int            **smaps;
   hypre_StructStencil ***sstencils;
   hypre_StructMatrix  ***smatrices;
   HYPRE_Int            **symmetric;
   HYPRE_Int            **num_centries;
   HYPRE_Int           ***centries;

   hypre_StructStencil   *sstencil;
   HYPRE_Int             *vars;
   hypre_Index           *sstencil_shape;
   HYPRE_Int              sstencil_size;
   HYPRE_Int             *new_sizes;
   hypre_Index          **new_shapes;
   HYPRE_Int              size;
   hypre_StructGrid      *sgrid;

   HYPRE_Int              vi, vj;
   HYPRE_Int              i, j, k;

   pmatrix = hypre_TAlloc(hypre_SStructPMatrix, 1);

   hypre_SStructPMatrixComm(pmatrix)     = comm;
   hypre_SStructPMatrixPGrid(pmatrix)    = pgrid;
   hypre_SStructPMatrixStencils(pmatrix) = stencils;
   hypre_SStructPMatrixNVars(pmatrix)    = nvars;

   /* create sstencils */
   smaps      = hypre_TAlloc(HYPRE_Int *, nvars);
   sstencils  = hypre_TAlloc(hypre_StructStencil **, nvars);
   new_sizes  = hypre_TAlloc(HYPRE_Int, nvars);
   new_shapes = hypre_TAlloc(hypre_Index *, nvars);
   size = 0;
   for (vi = 0; vi < nvars; vi++)
   {
      sstencils[vi] = hypre_TAlloc(hypre_StructStencil *, nvars);
      for (vj = 0; vj < nvars; vj++)
      {
         sstencils[vi][vj] = NULL;
         new_sizes[vj] = 0;
      }

      sstencil       = hypre_SStructStencilSStencil(stencils[vi]);
      vars           = hypre_SStructStencilVars(stencils[vi]);
      sstencil_shape = hypre_StructStencilShape(sstencil);
      sstencil_size  = hypre_StructStencilSize(sstencil);

      smaps[vi] = hypre_TAlloc(HYPRE_Int, sstencil_size);
      for (i = 0; i < sstencil_size; i++)
      {
         j = vars[i];
         new_sizes[j]++;
      }
      for (vj = 0; vj < nvars; vj++)
      {
         if (new_sizes[vj])
         {
            new_shapes[vj] = hypre_TAlloc(hypre_Index, new_sizes[vj]);
            new_sizes[vj] = 0;
         }
      }
      for (i = 0; i < sstencil_size; i++)
      {
         j = vars[i];
         k = new_sizes[j];
         hypre_CopyIndex(sstencil_shape[i], new_shapes[j][k]);
         smaps[vi][i] = k;
         new_sizes[j]++;
      }
      for (vj = 0; vj < nvars; vj++)
      {
         if (new_sizes[vj])
         {
            sstencils[vi][vj] =
               hypre_StructStencilCreate(ndim, new_sizes[vj], new_shapes[vj]);
         }
         size = hypre_max(size, new_sizes[vj]);
      }
   }
   hypre_SStructPMatrixSMaps(pmatrix)     = smaps;
   hypre_SStructPMatrixSStencils(pmatrix) = sstencils;
   hypre_TFree(new_sizes);
   hypre_TFree(new_shapes);

   /* create smatrices */
   smatrices = hypre_TAlloc(hypre_StructMatrix **, nvars);
   for (vi = 0; vi < nvars; vi++)
   {
      smatrices[vi] = hypre_TAlloc(hypre_StructMatrix *, nvars);
      for (vj = 0; vj < nvars; vj++)
      {
         smatrices[vi][vj] = NULL;
         if (sstencils[vi][vj] != NULL)
         {
            sgrid = hypre_SStructPGridSGrid(pgrid, vi);
            smatrices[vi][vj] =
               hypre_StructMatrixCreate(comm, sgrid, sstencils[vi][vj]);
         }
      }
   }
   hypre_SStructPMatrixSMatrices(pmatrix) = smatrices;

   /* create domain and range grid strides */
   hypre_SetIndex(hypre_SStructPMatrixDomainStride(pmatrix), 1);
   hypre_SetIndex(hypre_SStructPMatrixRangeStride(pmatrix), 1);

   /* create arrays */
   symmetric     = hypre_TAlloc(HYPRE_Int *, nvars);
   num_centries  = hypre_TAlloc(HYPRE_Int *, nvars);
   centries      = hypre_TAlloc(HYPRE_Int **, nvars);
   for (vi = 0; vi < nvars; vi++)
   {
      symmetric[vi]    = hypre_TAlloc(HYPRE_Int, nvars);
      num_centries[vi] = hypre_TAlloc(HYPRE_Int, nvars);
      centries[vi]     = hypre_TAlloc(HYPRE_Int *, nvars);
      for (vj = 0; vj < nvars; vj++)
      {
         symmetric[vi][vj]    = 0;
         num_centries[vi][vj] = 0;
      }
   }

   hypre_SStructPMatrixSymmetric(pmatrix) = symmetric;
   hypre_SStructPMatrixNumCEntries(pmatrix) = num_centries;
   hypre_SStructPMatrixCEntries(pmatrix) = centries;
   hypre_SStructPMatrixSEntriesSize(pmatrix) = size;
   hypre_SStructPMatrixSEntries(pmatrix) = hypre_TAlloc(HYPRE_Int, size);
   hypre_SStructPMatrixRefCount(pmatrix) = 1;

   *pmatrix_ptr = pmatrix;

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructPMatrixDestroy( hypre_SStructPMatrix *pmatrix )
{
   hypre_SStructStencil  **stencils;
   HYPRE_Int               nvars;
   HYPRE_Int             **smaps;
   hypre_StructStencil  ***sstencils;
   hypre_StructMatrix   ***smatrices;
   HYPRE_Int             **symmetric;
   HYPRE_Int             **num_centries;
   HYPRE_Int            ***centries;
   HYPRE_Int              *sentries;

   HYPRE_Int               vi, vj;

   if (pmatrix)
   {
      hypre_SStructPMatrixRefCount(pmatrix) --;
      if (hypre_SStructPMatrixRefCount(pmatrix) == 0)
      {
         stencils     = hypre_SStructPMatrixStencils(pmatrix);
         nvars        = hypre_SStructPMatrixNVars(pmatrix);
         smaps        = hypre_SStructPMatrixSMaps(pmatrix);
         sstencils    = hypre_SStructPMatrixSStencils(pmatrix);
         smatrices    = hypre_SStructPMatrixSMatrices(pmatrix);
         symmetric    = hypre_SStructPMatrixSymmetric(pmatrix);
         num_centries = hypre_SStructPMatrixNumCEntries(pmatrix);
         centries     = hypre_SStructPMatrixCEntries(pmatrix);
         sentries     = hypre_SStructPMatrixSEntries(pmatrix);

         for (vi = 0; vi < nvars; vi++)
         {
            HYPRE_SStructStencilDestroy(stencils[vi]);
            hypre_TFree(smaps[vi]);
            for (vj = 0; vj < nvars; vj++)
            {
               hypre_StructStencilDestroy(sstencils[vi][vj]);
               hypre_StructMatrixDestroy(smatrices[vi][vj]);
               hypre_TFree(centries[vi][vj]);
            }
            hypre_TFree(sstencils[vi]);
            hypre_TFree(smatrices[vi]);
            hypre_TFree(symmetric[vi]);
            hypre_TFree(num_centries[vi]);
            hypre_TFree(centries[vi]);
         }
         hypre_TFree(stencils);
         hypre_TFree(smaps);
         hypre_TFree(sstencils);
         hypre_TFree(smatrices);
         hypre_TFree(symmetric);
         hypre_TFree(num_centries);
         hypre_TFree(centries);
         hypre_TFree(sentries);
         hypre_TFree(pmatrix);
      }
   }

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/
HYPRE_Int
hypre_SStructPMatrixInitialize( hypre_SStructPMatrix *pmatrix )
{
   HYPRE_Int             nvars        = hypre_SStructPMatrixNVars(pmatrix);
   HYPRE_Int           **symmetric    = hypre_SStructPMatrixSymmetric(pmatrix);
   HYPRE_Int           **num_centries = hypre_SStructPMatrixNumCEntries(pmatrix);
   HYPRE_Int          ***centries     = hypre_SStructPMatrixCEntries(pmatrix);
   hypre_IndexRef        dom_stride   = hypre_SStructPMatrixDomainStride(pmatrix);
   hypre_IndexRef        ran_stride   = hypre_SStructPMatrixRangeStride(pmatrix);
   HYPRE_Int             num_ghost[2*HYPRE_MAXDIM];
   hypre_StructMatrix   *smatrix;
   HYPRE_Int             vi, vj, d, ndim;

   ndim = hypre_SStructPMatrixNDim(pmatrix);
   /* RDF: Why are the ghosts being reset to one? Maybe it needs to be at least
    * one to set shared coefficients correctly, but not exactly one? */
   for (d = 0; d < ndim; d++)
   {
      num_ghost[2*d] = num_ghost[2*d+1] = 1;
   }
   for (d = ndim; d < ndim; d++)
   {
      num_ghost[2*d] = num_ghost[2*d+1] = 0;
   }
   for (vi = 0; vi < nvars; vi++)
   {
      for (vj = 0; vj < nvars; vj++)
      {
         smatrix = hypre_SStructPMatrixSMatrix(pmatrix, vi, vj);
         if (smatrix != NULL)
         {
            HYPRE_StructMatrixSetDomainStride(smatrix, dom_stride);
            HYPRE_StructMatrixSetRangeStride(smatrix, ran_stride);
            HYPRE_StructMatrixSetConstantEntries(smatrix,
                                                 num_centries[vi][vj],
                                                 centries[vi][vj]);
            HYPRE_StructMatrixSetSymmetric(smatrix, symmetric[vi][vj]);
            HYPRE_StructMatrixSetNumGhost(smatrix, num_ghost);
            hypre_StructMatrixInitialize(smatrix);
            /* needed to get AddTo accumulation correct between processors */
            hypre_StructMatrixClearGhostValues(smatrix);
         }
      }
   }

   hypre_SStructPMatrixAccumulated(pmatrix) = 0;

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * (action > 0): add-to values
 * (action = 0): set values
 * (action < 0): get values
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructPMatrixSetValues( hypre_SStructPMatrix *pmatrix,
                               hypre_Index           index,
                               HYPRE_Int             var,
                               HYPRE_Int             nentries,
                               HYPRE_Int            *entries,
                               HYPRE_Complex        *values,
                               HYPRE_Int             action )
{
   hypre_SStructStencil *stencil = hypre_SStructPMatrixStencil(pmatrix, var);
   HYPRE_Int            *smap    = hypre_SStructPMatrixSMap(pmatrix, var);
   HYPRE_Int            *vars    = hypre_SStructStencilVars(stencil);
   hypre_StructMatrix   *smatrix;
   hypre_BoxArray       *grid_boxes;
   hypre_Box            *box, *grow_box;
   HYPRE_Int            *sentries;
   HYPRE_Int             i;

   smatrix = hypre_SStructPMatrixSMatrix(pmatrix, var, vars[entries[0]]);

   sentries = hypre_SStructPMatrixSEntries(pmatrix);
   for (i = 0; i < nentries; i++)
   {
      sentries[i] = smap[entries[i]];
   }

   /* set values inside the grid */
   hypre_StructMatrixSetValues(smatrix, index, nentries, sentries, values,
                               action, -1, 0);

   /* set (AddTo/Get) or clear (Set) values outside the grid in ghost zones */
   if (action != 0)
   {
      /* AddTo/Get */
      hypre_SStructPGrid  *pgrid = hypre_SStructPMatrixPGrid(pmatrix);
      hypre_Index          varoffset;
      HYPRE_Int            done = 0;

      grid_boxes = hypre_StructGridBoxes(hypre_StructMatrixGrid(smatrix));

      hypre_ForBoxI(i, grid_boxes)
      {
         box = hypre_BoxArrayBox(grid_boxes, i);
         if (hypre_IndexInBox(index, box))
         {
            done = 1;
            break;
         }
      }

      if (!done)
      {
         grow_box = hypre_BoxCreate(hypre_BoxArrayNDim(grid_boxes));
         hypre_SStructVariableGetOffset(hypre_SStructPGridVarType(pgrid, var),
                                        hypre_SStructPGridNDim(pgrid), varoffset);
         hypre_ForBoxI(i, grid_boxes)
         {
            box = hypre_BoxArrayBox(grid_boxes, i);
            hypre_CopyBox(box, grow_box);
            hypre_BoxGrowByIndex(grow_box, varoffset);
            if (hypre_IndexInBox(index, grow_box))
            {
               hypre_StructMatrixSetValues(smatrix, index, nentries, sentries,
                                           values, action, i, 1);
               break;
            }
         }
         hypre_BoxDestroy(grow_box);
      }
   }
   else
   {
      /* Set */
      grid_boxes = hypre_StructGridBoxes(hypre_StructMatrixGrid(smatrix));

      hypre_ForBoxI(i, grid_boxes)
      {
         box = hypre_BoxArrayBox(grid_boxes, i);
         if (!hypre_IndexInBox(index, box))
         {
            hypre_StructMatrixClearValues(smatrix, index, nentries, sentries, i, 1);
         }
      }
   }

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * (action > 0): add-to values
 * (action = 0): set values
 * (action < 0): get values
 * (action =-2): get values and zero out
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructPMatrixSetBoxValues( hypre_SStructPMatrix *pmatrix,
                                  hypre_Index           ilower,
                                  hypre_Index           iupper,
                                  HYPRE_Int             var,
                                  HYPRE_Int             nentries,
                                  HYPRE_Int            *entries,
                                  HYPRE_Complex        *values,
                                  HYPRE_Int             action )
{
   HYPRE_Int             ndim    = hypre_SStructPMatrixNDim(pmatrix);
   hypre_SStructStencil *stencil = hypre_SStructPMatrixStencil(pmatrix, var);
   HYPRE_Int            *smap    = hypre_SStructPMatrixSMap(pmatrix, var);
   HYPRE_Int            *vars    = hypre_SStructStencilVars(stencil);
   hypre_StructMatrix   *smatrix;
   hypre_BoxArray       *grid_boxes;
   hypre_Box            *box;
   hypre_Box            *value_box;
   HYPRE_Int            *sentries;
   HYPRE_Int             i, j;

   smatrix = hypre_SStructPMatrixSMatrix(pmatrix, var, vars[entries[0]]);

   box = hypre_BoxCreate(hypre_StructMatrixNDim(smatrix));
   hypre_CopyIndex(ilower, hypre_BoxIMin(box));
   hypre_CopyIndex(iupper, hypre_BoxIMax(box));
   value_box = box;

   sentries = hypre_SStructPMatrixSEntries(pmatrix);
   for (i = 0; i < nentries; i++)
   {
      sentries[i] = smap[entries[i]];
   }

   /* set values inside the grid */
   hypre_StructMatrixSetBoxValues(smatrix, box, value_box, nentries, sentries,
                                  values, action, -1, 0);

   /* set (AddTo/Get) or clear (Set) values outside the grid in ghost zones */
   if (action != 0)
   {
      /* AddTo/Get */
      hypre_SStructPGrid  *pgrid = hypre_SStructPMatrixPGrid(pmatrix);
      hypre_Index          varoffset;
      hypre_BoxArray      *left_boxes, *done_boxes, *temp_boxes;
      hypre_Box           *left_box, *done_box, *int_box;

      hypre_SStructVariableGetOffset(hypre_SStructPGridVarType(pgrid, var),
                                     hypre_SStructPGridNDim(pgrid), varoffset);
      grid_boxes = hypre_StructGridBoxes(hypre_StructMatrixGrid(smatrix));

      left_boxes = hypre_BoxArrayCreate(1, ndim);
      done_boxes = hypre_BoxArrayCreate(2, ndim);
      temp_boxes = hypre_BoxArrayCreate(0, ndim);

      /* done_box always points to the first box in done_boxes */
      done_box = hypre_BoxArrayBox(done_boxes, 0);
      /* int_box always points to the second box in done_boxes */
      int_box = hypre_BoxArrayBox(done_boxes, 1);

      hypre_CopyBox(box, hypre_BoxArrayBox(left_boxes, 0));
      hypre_BoxArraySetSize(left_boxes, 1);
      hypre_SubtractBoxArrays(left_boxes, grid_boxes, temp_boxes);

      hypre_BoxArraySetSize(done_boxes, 0);
      hypre_ForBoxI(i, grid_boxes)
      {
         hypre_SubtractBoxArrays(left_boxes, done_boxes, temp_boxes);
         hypre_BoxArraySetSize(done_boxes, 1);
         hypre_CopyBox(hypre_BoxArrayBox(grid_boxes, i), done_box);
         hypre_BoxGrowByIndex(done_box, varoffset);
         hypre_ForBoxI(j, left_boxes)
         {
            left_box = hypre_BoxArrayBox(left_boxes, j);
            hypre_IntersectBoxes(left_box, done_box, int_box);
            hypre_StructMatrixSetBoxValues(smatrix, int_box, value_box,
                                           nentries, sentries,
                                           values, action, i, 1);
         }
      }

      hypre_BoxArrayDestroy(left_boxes);
      hypre_BoxArrayDestroy(done_boxes);
      hypre_BoxArrayDestroy(temp_boxes);
   }
   else
   {
      /* Set */
      hypre_BoxArray  *diff_boxes;
      hypre_Box       *grid_box, *diff_box;

      grid_boxes = hypre_StructGridBoxes(hypre_StructMatrixGrid(smatrix));
      diff_boxes = hypre_BoxArrayCreate(0, ndim);

      hypre_ForBoxI(i, grid_boxes)
      {
         grid_box = hypre_BoxArrayBox(grid_boxes, i);
         hypre_BoxArraySetSize(diff_boxes, 0);
         hypre_SubtractBoxes(box, grid_box, diff_boxes);

         hypre_ForBoxI(j, diff_boxes)
         {
            diff_box = hypre_BoxArrayBox(diff_boxes, j);
            hypre_StructMatrixClearBoxValues(smatrix, diff_box, nentries, sentries,
                                             i, 1);
         }
      }
      hypre_BoxArrayDestroy(diff_boxes);
   }

   hypre_BoxDestroy(box);

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructPMatrixAccumulate( hypre_SStructPMatrix *pmatrix )
{
   hypre_SStructPGrid    *pgrid    = hypre_SStructPMatrixPGrid(pmatrix);
   HYPRE_Int              nvars    = hypre_SStructPMatrixNVars(pmatrix);
   HYPRE_Int              ndim     = hypre_SStructPGridNDim(pgrid);
   HYPRE_SStructVariable *vartypes = hypre_SStructPGridVarTypes(pgrid);

   hypre_StructMatrix    *smatrix;
   hypre_Index            varoffset;
   HYPRE_Int              num_ghost[2*HYPRE_MAXDIM];
   hypre_StructGrid      *sgrid;
   HYPRE_Int              vi, vj, d;

   hypre_CommInfo        *comm_info;
   hypre_CommPkg         *comm_pkg;
   hypre_CommHandle      *comm_handle;
   HYPRE_Complex         *data;

   /* if values already accumulated, just return */
   if (hypre_SStructPMatrixAccumulated(pmatrix))
   {
      return hypre_error_flag;
   }

   for (d = ndim; d < ndim; d++)
   {
      num_ghost[2*d] = num_ghost[2*d+1] = 0;
   }
   for (vi = 0; vi < nvars; vi++)
   {
      for (vj = 0; vj < nvars; vj++)
      {
         smatrix = hypre_SStructPMatrixSMatrix(pmatrix, vi, vj);
         if (smatrix != NULL)
         {
            sgrid = hypre_StructMatrixGrid(smatrix);
            /* assumes vi and vj vartypes are the same */
            hypre_SStructVariableGetOffset(vartypes[vi], ndim, varoffset);
            for (d = 0; d < ndim; d++)
            {
               num_ghost[2*d]   = num_ghost[2*d+1] = hypre_IndexD(varoffset, d);
            }

            /* accumulate values from AddTo */
            hypre_CreateCommInfoFromNumGhost(sgrid, num_ghost, &comm_info);
            hypre_CommPkgCreate(comm_info,
                                hypre_StructMatrixDataSpace(smatrix),
                                hypre_StructMatrixDataSpace(smatrix),
                                hypre_StructMatrixNumValues(smatrix), NULL, 1,
                                hypre_StructMatrixComm(smatrix),
                                &comm_pkg);
            data = hypre_StructMatrixVData(smatrix);
            hypre_InitializeCommunication(comm_pkg, &data, &data, 1, 0,
                                          &comm_handle);
            hypre_FinalizeCommunication(comm_handle);

            hypre_CommInfoDestroy(comm_info);
            hypre_CommPkgDestroy(comm_pkg);
         }
      }
   }

   hypre_SStructPMatrixAccumulated(pmatrix) = 1;

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructPMatrixAssemble( hypre_SStructPMatrix *pmatrix )
{
   HYPRE_Int              nvars    = hypre_SStructPMatrixNVars(pmatrix);
   hypre_StructMatrix    *smatrix;
   HYPRE_Int              vi, vj;

   hypre_SStructPMatrixAccumulate(pmatrix);

   for (vi = 0; vi < nvars; vi++)
   {
      for (vj = 0; vj < nvars; vj++)
      {
         smatrix = hypre_SStructPMatrixSMatrix(pmatrix, vi, vj);
         if (smatrix != NULL)
         {
            hypre_StructMatrixClearGhostValues(smatrix);
            hypre_StructMatrixAssemble(smatrix);
         }
      }
   }

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * TODO: Deprecate this function. var == -1 or to_var == -1 are never used.
 *       These cases are used only in HYPRE_SStructMatrixSetSymmetric.
 *--------------------------------------------------------------------------*/
#if 0
HYPRE_Int
hypre_SStructPMatrixSetSymmetric( hypre_SStructPMatrix *pmatrix,
                                  HYPRE_Int             var,
                                  HYPRE_Int             to_var,
                                  HYPRE_Int             symmetric )
{
   HYPRE_Int **pmsymmetric = hypre_SStructPMatrixSymmetric(pmatrix);

   HYPRE_Int vstart = var;
   HYPRE_Int vsize  = 1;
   HYPRE_Int tstart = to_var;
   HYPRE_Int tsize  = 1;
   HYPRE_Int v, t;

   if (var == -1)
   {
      vstart = 0;
      vsize  = hypre_SStructPMatrixNVars(pmatrix);
   }
   if (to_var == -1)
   {
      tstart = 0;
      tsize  = hypre_SStructPMatrixNVars(pmatrix);
   }

   for (v = vstart; v < vsize; v++)
   {
      for (t = tstart; t < tsize; t++)
      {
         pmsymmetric[v][t] = symmetric;
      }
   }

   return hypre_error_flag;
}

#else
/*--------------------------------------------------------------------------
 * NOTE: Should we have an accessor macro for doing this job?
 *       How would we call it?
 *--------------------------------------------------------------------------*/
HYPRE_Int
hypre_SStructPMatrixSetSymmetric( hypre_SStructPMatrix *pmatrix,
                                  HYPRE_Int             var,
                                  HYPRE_Int             to_var,
                                  HYPRE_Int             symmetric )
{
   HYPRE_Int **pmsymmetric = hypre_SStructPMatrixSymmetric(pmatrix);

   pmsymmetric[var][to_var] = symmetric;

   return hypre_error_flag;
}
#endif

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/
HYPRE_Int
hypre_SStructPMatrixSetCEntries( hypre_SStructPMatrix *pmatrix,
                                 HYPRE_Int             var,
                                 HYPRE_Int             to_var,
                                 HYPRE_Int             num_centries,
                                 HYPRE_Int            *centries )
{
   HYPRE_Int   **pmnum_centries = hypre_SStructPMatrixNumCEntries(pmatrix);
   HYPRE_Int  ***pmcentries     = hypre_SStructPMatrixCEntries(pmatrix);
   HYPRE_Int     i;

   pmnum_centries[var][to_var] = num_centries;
   pmcentries[var][to_var]     = hypre_CTAlloc(HYPRE_Int, num_centries);
   for (i = 0; i < num_centries; i++)
   {
      pmcentries[var][to_var][i] = centries[i];
   }

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * NOTE: Should we have an accessor macro for doing this job?
 *       How would we call it?
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructPMatrixSetDomainStride( hypre_SStructPMatrix *pmatrix,
                                     hypre_Index           dom_stride )
{
   hypre_CopyIndex(dom_stride, hypre_SStructPMatrixDomainStride(pmatrix));

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * NOTE: Should we have an accessor macro for doing this job?
 *       How would we call it?
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructPMatrixSetRangeStride( hypre_SStructPMatrix *pmatrix,
                                    hypre_Index           ran_stride )
{
   hypre_CopyIndex(ran_stride, hypre_SStructPMatrixRangeStride(pmatrix));

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructPMatrixPrint( const char           *filename,
                           hypre_SStructPMatrix *pmatrix,
                           HYPRE_Int             all )
{
   HYPRE_Int           nvars = hypre_SStructPMatrixNVars(pmatrix);
   hypre_StructMatrix *smatrix;
   HYPRE_Int           vi, vj;
   char                new_filename[255];

   for (vi = 0; vi < nvars; vi++)
   {
      for (vj = 0; vj < nvars; vj++)
      {
         smatrix = hypre_SStructPMatrixSMatrix(pmatrix, vi, vj);
         if (smatrix != NULL)
         {
            hypre_sprintf(new_filename, "%s.v%1d%1d", filename, vi, vj);
            hypre_StructMatrixPrint(new_filename, smatrix, all);
         }
      }
   }

   return hypre_error_flag;
}

/*==========================================================================
 * SStructUMatrix routines
 *==========================================================================*/

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructUMatrixInitialize( hypre_SStructMatrix *matrix )
{
   HYPRE_Int               ndim          = hypre_SStructMatrixNDim(matrix);
   HYPRE_IJMatrix          ijmatrix      = hypre_SStructMatrixIJMatrix(matrix);
   HYPRE_Int               matrix_type   = hypre_SStructMatrixObjectType(matrix);
   hypre_SStructGraph     *graph         = hypre_SStructMatrixGraph(matrix);
   hypre_SStructStencil ***stencils      = hypre_SStructGraphStencils(graph);
   HYPRE_Int               nUventries    = hypre_SStructGraphNUVEntries(graph);
   HYPRE_Int              *iUventries    = hypre_SStructGraphIUVEntries(graph);
   hypre_SStructUVEntry  **Uventries     = hypre_SStructGraphUVEntries(graph);
   HYPRE_Int              *active_pids   = hypre_SStructGraphActivePartIDs(graph);
   HYPRE_Int               active_nparts = hypre_SStructGraphActiveNParts(graph);

   hypre_SStructGrid      *ran_grid      = hypre_SStructGraphRanGrid(graph);
   HYPRE_Int              *ran_pids      = hypre_SStructGridPartIDs(ran_grid);
   HYPRE_Int               ran_nparts    = hypre_SStructGridNParts(ran_grid);
   HYPRE_Int             **nvneighbors   = hypre_SStructGridNVNeighbors(ran_grid);

   hypre_SStructPGrid     *pgrid;
   hypre_StructGrid       *sgrid;
   hypre_SStructStencil   *stencil;
   hypre_BoxArray         *boxes;
   hypre_Box              *box;
   hypre_Box              *ghost_box;
   hypre_IndexRef          start;
   hypre_Index             loop_size, stride;

   HYPRE_Int              *split;
   HYPRE_Int               nvars;
   HYPRE_Int               nrows, rowstart, nnzs;
   HYPRE_Int               part, ran_part, var, entry, b, m, mi;
   HYPRE_Int              *row_sizes;
   HYPRE_Int               max_size;

   HYPRE_IJMatrixSetObjectType(ijmatrix, HYPRE_PARCSR);
   if (matrix_type == HYPRE_SSTRUCT || matrix_type == HYPRE_STRUCT)
   {
      rowstart = hypre_SStructGridGhstartRank(ran_grid);
      nrows = hypre_SStructGridGhlocalSize(ran_grid);
   }
   else /* matrix_type == HYPRE_PARCSR */
   {
      rowstart = hypre_SStructGridStartRank(ran_grid);
      nrows = hypre_SStructGridLocalSize(ran_grid);
   }

   /* set row sizes */
   max_size  = m = 0;
   ghost_box = hypre_BoxCreate(ndim);
   row_sizes = hypre_CTAlloc(HYPRE_Int, nrows);
   hypre_SetIndex(stride, 1);
   for (ran_part = 0; ran_part < ran_nparts; ran_part++)
   {
      part  = hypre_BinarySearch(active_pids, ran_pids[ran_part], active_nparts);
      pgrid = hypre_SStructGridPGrid(ran_grid, ran_part);
      nvars = hypre_SStructPGridNVars(pgrid);

      if (part > -1)
      {
         /* This part is active in the range grid */
         for (var = 0; var < nvars; var++)
         {
            sgrid = hypre_SStructPGridSGrid(pgrid, var);

            stencil = stencils[part][var];
            split = hypre_SStructMatrixSplit(matrix, part, var);
            nnzs = 0;
            for (entry = 0; entry < hypre_SStructStencilSize(stencil); entry++)
            {
               if (split[entry] == -1)
               {
                  nnzs++;
               }
            }
#if 0
            /* TODO: For now, assume stencil is full/complete */
            if (hypre_SStructMatrixSymmetric(matrix))
            {
               nnzs = 2*nnzs - 1;
            }
#endif
            boxes = hypre_StructGridBoxes(sgrid);
            hypre_ForBoxI(b, boxes)
            {
               box = hypre_BoxArrayBox(boxes, b);
               hypre_CopyBox(box, ghost_box);
               if (matrix_type == HYPRE_SSTRUCT || matrix_type == HYPRE_STRUCT)
               {
                  hypre_BoxGrowByArray(ghost_box, hypre_StructGridNumGhost(sgrid));
               }

               start = hypre_BoxIMin(box);
               hypre_BoxGetSize(box, loop_size);
               hypre_BoxLoop1Begin(ndim, loop_size, ghost_box, start, stride, mi);
#ifdef HYPRE_USING_OPENMP
#pragma omp parallel for private(HYPRE_BOX_PRIVATE,mi) HYPRE_SMP_SCHEDULE
#endif
               hypre_BoxLoop1For(mi)
               {
                  row_sizes[m + mi] = nnzs;
               }
               hypre_BoxLoop1End(mi);

               m += hypre_BoxVolume(ghost_box);
            }

            max_size = hypre_max(max_size, nnzs);
            if (nvneighbors[part][var])
            {
               max_size = hypre_max(max_size, hypre_SStructStencilSize(stencil));
            }
         } /* loop on variables */
      }
      else
      {
         /* This part is not active in the range grid */
         for (var = 0; var < nvars; var++)
         {
            sgrid = hypre_SStructPGridSGrid(pgrid, var);
            boxes = hypre_StructGridBoxes(sgrid);
            hypre_ForBoxI(b, boxes)
            {
               box = hypre_BoxArrayBox(boxes, b);
               hypre_CopyBox(box, ghost_box);
               if (matrix_type == HYPRE_SSTRUCT || matrix_type == HYPRE_STRUCT)
               {
                  hypre_BoxGrowByArray(ghost_box, hypre_StructGridNumGhost(sgrid));
               }

               start = hypre_BoxIMin(box);
               hypre_BoxGetSize(box, loop_size);
               hypre_BoxLoop1Begin(ndim, loop_size, ghost_box, start, stride, mi);
#ifdef HYPRE_USING_OPENMP
#pragma omp parallel for private(HYPRE_BOX_PRIVATE,mi) HYPRE_SMP_SCHEDULE
#endif
               hypre_BoxLoop1For(mi)
               {
                  row_sizes[m + mi] = 1;
               }
               hypre_BoxLoop1End(mi);

               m += hypre_BoxVolume(ghost_box);
            }
         } /* loop on variables */
      } /* if (part > -1) */
   } /* loop on parts */
   hypre_BoxDestroy(ghost_box);

   /* GEC0902 essentially for each UVentry we figure out how many extra columns
    * we need to add to the rowsizes                                   */

   /* RDF: THREAD? */
   for (entry = 0; entry < nUventries; entry++)
   {
      mi = iUventries[entry];
      m = hypre_SStructUVEntryRank(Uventries[mi]) - rowstart;
      if ((m > -1) && (m < nrows))
      {
         row_sizes[m] += hypre_SStructUVEntryNUEntries(Uventries[mi]);
         max_size = hypre_max(max_size, row_sizes[m]);
      }
   }

   /* ZTODO: Update row_sizes based on neighbor off-part couplings */
   HYPRE_IJMatrixSetRowSizes(ijmatrix, (const HYPRE_Int *) row_sizes);

   hypre_TFree(row_sizes);
   hypre_SStructMatrixTmpColCoords(matrix) = hypre_CTAlloc(HYPRE_Int, max_size);
   hypre_SStructMatrixTmpCoeffs(matrix) = hypre_CTAlloc(HYPRE_Complex, max_size);

   HYPRE_IJMatrixInitialize(ijmatrix);

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * (action > 0): add-to values
 * (action = 0): set values
 * (action < 0): get values
 *
 * 9/09 - AB: modified to use the box manager - here we need to check the
 *            neighbor box manager also
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructUMatrixSetValues( hypre_SStructMatrix *matrix,
                               HYPRE_Int            part,
                               hypre_Index          index,
                               HYPRE_Int            var,
                               HYPRE_Int            nentries,
                               HYPRE_Int           *entries,
                               HYPRE_Complex       *values,
                               HYPRE_Int            action )
{
   HYPRE_Int                ndim     = hypre_SStructMatrixNDim(matrix);
   HYPRE_IJMatrix           ijmatrix = hypre_SStructMatrixIJMatrix(matrix);
   hypre_SStructGraph      *graph    = hypre_SStructMatrixGraph(matrix);
   hypre_SStructGrid       *dom_grid = hypre_SStructGraphDomGrid(graph);
   hypre_SStructGrid       *ran_grid = hypre_SStructGraphRanGrid(graph);
   hypre_SStructStencil    *stencil  = hypre_SStructGraphStencil(graph, part, var);
   HYPRE_Int               *pmaps    = hypre_SStructGraphActivePMaps(graph);
   HYPRE_Int               *vars     = hypre_SStructStencilVars(stencil);
   hypre_Index             *shape    = hypre_SStructStencilShape(stencil);
   HYPRE_Int                size     = hypre_SStructStencilSize(stencil);
   hypre_IndexRef           offset;
   hypre_Index              to_index;
   hypre_SStructUVEntry    *Uventry;
   hypre_BoxManEntry       *boxman_entry;
   hypre_SStructBoxManInfo *entry_info;
   HYPRE_Int                row_coord;
   HYPRE_Int               *col_coords;
   HYPRE_Int                ncoeffs;
   HYPRE_Complex           *coeffs;
   HYPRE_Int                i, entry, Uverank;
   HYPRE_Int                matrix_type = hypre_SStructMatrixObjectType(matrix);

   hypre_SStructGridFindBoxManEntry(ran_grid, pmaps[part], index, var, &boxman_entry);

   /* if not local, check neighbors */
   if (boxman_entry == NULL)
   {
      hypre_SStructGridFindNborBoxManEntry(ran_grid, pmaps[part], index, var,
                                           &boxman_entry);
   }

   if (boxman_entry == NULL)
   {
      hypre_error_in_arg(1);
      hypre_error_in_arg(2);
      hypre_error_in_arg(3);
      return hypre_error_flag;
   }
   else
   {
      hypre_BoxManEntryGetInfo(boxman_entry, (void **) &entry_info);
   }

   hypre_SStructBoxManEntryGetGlobalRank(boxman_entry, index,
                                         &row_coord, matrix_type);

   col_coords = hypre_SStructMatrixTmpColCoords(matrix);
   coeffs     = hypre_SStructMatrixTmpCoeffs(matrix);
   ncoeffs = 0;
   for (i = 0; i < nentries; i++)
   {
      entry = entries[i];

      if (entry < size)
      {
         /* stencil entries */
         offset = shape[entry];
         hypre_AddIndexes(index, offset, ndim, to_index);

         hypre_SStructGridFindBoxManEntry(dom_grid, part, to_index, vars[entry],
                                          &boxman_entry);

         /* if not local, check neighbors */
         if (boxman_entry == NULL)
            hypre_SStructGridFindNborBoxManEntry(dom_grid, part, to_index,
                                                 vars[entry], &boxman_entry);

         if (boxman_entry != NULL)
         {
            hypre_SStructBoxManEntryGetGlobalRank(boxman_entry, to_index,
                                                  &col_coords[ncoeffs],matrix_type);

            coeffs[ncoeffs] = values[i];
            ncoeffs++;
         }
      }
      else
      {
         /* non-stencil entries */
         entry -= size;
         hypre_SStructGraphGetUVEntryRank(graph, part, var, index, &Uverank);

         if (Uverank > -1)
         {
            Uventry = hypre_SStructGraphUVEntry(graph, Uverank);
            col_coords[ncoeffs] = hypre_SStructUVEntryToRank(Uventry, entry);
            coeffs[ncoeffs] = values[i];
            ncoeffs++;
         }
      }
   }

   if (action > 0)
   {
      HYPRE_IJMatrixAddToValues(ijmatrix, 1, &ncoeffs, &row_coord,
                                (const HYPRE_Int *) col_coords,
                                (const HYPRE_Complex *) coeffs);
   }
   else if (action > -1)
   {
      HYPRE_IJMatrixSetValues(ijmatrix, 1, &ncoeffs, &row_coord,
                              (const HYPRE_Int *) col_coords,
                              (const HYPRE_Complex *) coeffs);
   }
   else
   {
      HYPRE_IJMatrixGetValues(ijmatrix, 1, &ncoeffs, &row_coord,
                              col_coords, values);
   }

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * Note: Entries must all be of type stencil or non-stencil, but not both.
 *
 * (action > 0): add-to values
 * (action = 0): set values
 * (action < 0): get values
 *
 * 9/09 - AB: modified to use the box manager- here we need to check the
 *            neighbor box manager also
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructUMatrixSetBoxValuesHelper( hypre_SStructMatrix *matrix,
                                        HYPRE_Int            part,
                                        hypre_Index          ilower,
                                        hypre_Index          iupper,
                                        HYPRE_Int            var,
                                        HYPRE_Int            nentries,
                                        HYPRE_Int           *entries,
                                        HYPRE_Complex       *values,
                                        HYPRE_Int            action,
                                        HYPRE_IJMatrix       ijmatrix )
{
   HYPRE_Int             ndim        = hypre_SStructMatrixNDim(matrix);
   HYPRE_Int             matrix_type = hypre_SStructMatrixObjectType(matrix);
   hypre_SStructPMatrix *pmatrix     = hypre_SStructMatrixPMatrix(matrix, part);
   hypre_SStructGraph   *graph       = hypre_SStructMatrixGraph(matrix);
   hypre_SStructGrid    *ran_grid    = hypre_SStructGraphRanGrid(graph);
   hypre_SStructGrid    *dom_grid    = hypre_SStructGraphDomGrid(graph);
   hypre_SStructStencil *stencil     = hypre_SStructGraphStencil(graph, part, var);
   HYPRE_Int            *pmaps       = hypre_SStructGraphActivePMaps(graph);
   HYPRE_Int            *vars        = hypre_SStructStencilVars(stencil);
   hypre_Index          *shape       = hypre_SStructStencilShape(stencil);
   HYPRE_Int             size        = hypre_SStructStencilSize(stencil);
   hypre_IndexRef        dom_stride  = hypre_SStructPMatrixDomainStride(pmatrix);
   hypre_IndexRef        offset;
   hypre_BoxManEntry   **boxman_entries;
   HYPRE_Int             nboxman_entries;
   hypre_BoxManEntry   **boxman_to_entries;
   HYPRE_Int             nboxman_to_entries;
   HYPRE_Int             nrows;
   HYPRE_Int            *ncols;
   HYPRE_Int            *rows;
   HYPRE_Int            *cols;
   HYPRE_Complex        *ijvalues;
   hypre_Box            *box, *vbox;
   hypre_Box            *to_box;
   hypre_Box            *map_box;
   hypre_Box            *int_box;
   hypre_Box            *map_vbox;
   hypre_Index           index, unit_stride, loop_size;
   hypre_IndexRef        start;
   hypre_Index           origin;
   hypre_Index           rs, cs;
   HYPRE_Int             row_base, col_base;
   HYPRE_Int             d, ei, entry, ii, jj, i, mi, vi;

   box  = hypre_BoxCreate(ndim);
   vbox = hypre_BoxCreate(ndim);

   hypre_BoxSetExtents(vbox, ilower, iupper);

   /*------------------------------------------
    * all stencil entries
    *------------------------------------------*/

   if (entries[0] < size)
   {
      to_box   = hypre_BoxCreate(ndim);
      map_box  = hypre_BoxCreate(ndim);
      int_box  = hypre_BoxCreate(ndim);
      map_vbox = hypre_BoxCreate(ndim);

      nrows    = hypre_BoxVolume(vbox)*nentries;
      ncols    = hypre_CTAlloc(HYPRE_Int, nrows);
#ifdef HYPRE_USING_OPENMP
#pragma omp parallel for private(i) HYPRE_SMP_SCHEDULE
#endif
      for (i = 0; i < nrows; i++)
      {
         ncols[i] = 1;
      }
      rows     = hypre_CTAlloc(HYPRE_Int, nrows);
      cols     = hypre_CTAlloc(HYPRE_Int, nrows);
      ijvalues = hypre_CTAlloc(HYPRE_Complex, nrows);

      hypre_SetIndex(unit_stride, 1);
      hypre_SetIndex(origin, 0);

      hypre_SStructGridIntersect(ran_grid, pmaps[part], var, vbox, -1,
                                 &boxman_entries, &nboxman_entries);

      for (ii = 0; ii < nboxman_entries; ii++)
      {
         hypre_SStructBoxManEntryGetStrides(boxman_entries[ii], rs, matrix_type);

         hypre_CopyBox(vbox, box);
         hypre_BoxManEntryGetExtents(boxman_entries[ii],
                                     hypre_BoxIMin(map_box),
                                     hypre_BoxIMax(map_box));
         hypre_IntersectBoxes(box, map_box, int_box);
         hypre_CopyBox(int_box, box);

         nrows = 0;
         for (ei = 0; ei < nentries; ei++)
         {
            entry = entries[ei];
            offset = shape[entry];

            hypre_CopyBox(box, to_box);
            hypre_BoxShiftPos(to_box, offset);
            hypre_CoarsenBox(to_box, NULL, dom_stride);

            hypre_SStructGridIntersect(dom_grid, part, vars[entry], to_box, -1,
                                       &boxman_to_entries, &nboxman_to_entries);

            for (jj = 0; jj < nboxman_to_entries; jj++)
            {
               hypre_SStructBoxManEntryGetStrides(boxman_to_entries[jj],
                                                  cs, matrix_type);

               hypre_BoxManEntryGetExtents(boxman_to_entries[jj],
                                           hypre_BoxIMin(map_box),
                                           hypre_BoxIMax(map_box));
               hypre_IntersectBoxes(to_box, map_box, int_box);

               hypre_CopyIndex(hypre_BoxIMin(int_box), index);
               hypre_SStructBoxManEntryGetGlobalRank(boxman_to_entries[jj],
                                                     index, &col_base, matrix_type);

               hypre_RefineBox(int_box, NULL, dom_stride);
               hypre_BoxShiftNeg(int_box, offset);

               hypre_CopyIndex(hypre_BoxIMin(int_box), index);
               hypre_SStructBoxManEntryGetGlobalRank(boxman_entries[ii],
                                                     index, &row_base, matrix_type);
               hypre_CopyBox(vbox, map_vbox);

               hypre_SStructMatrixMapDataBox(matrix, part, var, vars[entry], map_vbox);
               hypre_SStructMatrixMapDataBox(matrix, part, var, vars[entry], int_box);

               start = hypre_BoxIMin(int_box);
               hypre_BoxGetSize(int_box, loop_size);

               hypre_BoxLoop2Begin(ndim, loop_size,
                                   int_box,  start, unit_stride, mi,
                                   map_vbox, start, unit_stride, vi);
#ifdef HYPRE_USING_OPENMP
#pragma omp parallel for private(HYPRE_BOX_PRIVATE,mi,vi,index,d) HYPRE_SMP_SCHEDULE
#endif
               hypre_BoxLoop2For(mi, vi)
               {
                  hypre_BoxLoopGetIndex(index);
                  rows[nrows + mi] = row_base;
                  cols[nrows + mi] = col_base;
                  for (d = 0; d < ndim; d++)
                  {
                     rows[nrows + mi] += index[d]*rs[d]*dom_stride[d];
                     cols[nrows + mi] += index[d]*cs[d];
                  }
                  ijvalues[nrows + mi] = values[ei + vi*nentries]; // TODO: change access pattern?
                  /* hypre_printf("[%d, %d, %d]: (%d, %d): %f\n", nrows, mi, (nrows + mi), */
                  /*               rows[nrows + mi], cols[nrows + mi], ijvalues[nrows + mi]); */
               }
               hypre_BoxLoop2End(mi, vi);

               nrows += hypre_BoxVolume(int_box);

            } /* end loop through boxman to entries */

            hypre_TFree(boxman_to_entries);

         } /* end of ei nentries loop */

         /*------------------------------------------
          * set IJ values one stencil entry at a time
          *------------------------------------------*/

         if (action > 0)
         {
            HYPRE_IJMatrixAddToValues(ijmatrix, nrows, ncols,
                                      (const HYPRE_Int *) rows,
                                      (const HYPRE_Int *) cols,
                                      (const HYPRE_Complex *) ijvalues);
         }
         else if (action > -1)
         {
            HYPRE_IJMatrixSetValues(ijmatrix, nrows, ncols,
                                    (const HYPRE_Int *) rows,
                                    (const HYPRE_Int *) cols,
                                    (const HYPRE_Complex *) ijvalues);
         }
         else
         {
            HYPRE_IJMatrixGetValues(ijmatrix, nrows, ncols, rows, cols, values);
         }
      } /* end loop through boxman entries */

      hypre_TFree(boxman_entries);

      hypre_TFree(ncols);
      hypre_TFree(rows);
      hypre_TFree(cols);
      hypre_TFree(ijvalues);

      hypre_BoxDestroy(to_box);
      hypre_BoxDestroy(map_box);
      hypre_BoxDestroy(int_box);
      hypre_BoxDestroy(map_vbox);
   }

   /*------------------------------------------
    * non-stencil entries
    *------------------------------------------*/

   else
   {
      /* RDF: THREAD (Check safety on UMatrixSetValues call) */
      hypre_BoxGetSize(vbox, loop_size);
      hypre_BoxLoop0Begin(ndim, loop_size);
      hypre_BoxLoopSetOneBlock();
      hypre_BoxLoop0For()
      {
         hypre_BoxLoopGetIndex(index);
         hypre_SStructUMatrixSetValues(matrix, part, index, var,
                                       nentries, entries, values, action);
         values += nentries;
      }
      hypre_BoxLoop0End();
   }

   hypre_BoxDestroy(box);
   hypre_BoxDestroy(vbox);

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/


HYPRE_Int
hypre_SStructUMatrixSetBoxValues( hypre_SStructMatrix *matrix,
                                  HYPRE_Int            part,
                                  hypre_Index          ilower,
                                  hypre_Index          iupper,
                                  HYPRE_Int            var,
                                  HYPRE_Int            nentries,
                                  HYPRE_Int           *entries,
                                  HYPRE_Complex       *values,
                                  HYPRE_Int            action )
{
   HYPRE_IJMatrix  ijmatrix  = hypre_SStructMatrixIJMatrix(matrix);

   hypre_SStructUMatrixSetBoxValuesHelper(matrix, part, ilower, iupper,
                                          var, nentries, entries, values,
                                          action, ijmatrix);

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructUMatrixAssemble( hypre_SStructMatrix *matrix )
{
   HYPRE_IJMatrix ijmatrix = hypre_SStructMatrixIJMatrix(matrix);

   HYPRE_IJMatrixAssemble(ijmatrix);
   HYPRE_IJMatrixGetObject(
      ijmatrix, (void **) &hypre_SStructMatrixParCSRMatrix(matrix));

   return hypre_error_flag;
}

/*==========================================================================
 * SStructMatrix routines
 *==========================================================================*/

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/
HYPRE_Int
hypre_SStructMatrixMapDataBox( hypre_SStructMatrix  *matrix,
                               HYPRE_Int             part,
                               HYPRE_Int             vi,
                               HYPRE_Int             vj,
                               hypre_Box            *map_vbox )
{
   HYPRE_Int             matrix_type = hypre_SStructMatrixObjectType(matrix);
   hypre_SStructPMatrix *pmatrix;
   hypre_StructMatrix   *smatrix;

   if (matrix_type == HYPRE_SSTRUCT || matrix_type == HYPRE_STRUCT)
   {
      pmatrix = hypre_SStructMatrixPMatrix(matrix, part);
      smatrix = hypre_SStructPMatrixSMatrix(pmatrix, vi, vj);
      hypre_StructMatrixMapDataBox(smatrix, map_vbox);
   }

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructMatrixRef( hypre_SStructMatrix  *matrix,
                        hypre_SStructMatrix **matrix_ref )
{
   hypre_SStructMatrixRefCount(matrix) ++;
   *matrix_ref = matrix;

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructMatrixSplitEntries( hypre_SStructMatrix *matrix,
                                 HYPRE_Int            part,
                                 HYPRE_Int            var,
                                 HYPRE_Int            nentries,
                                 HYPRE_Int           *entries,
                                 HYPRE_Int           *nSentries_ptr,
                                 HYPRE_Int          **Sentries_ptr,
                                 HYPRE_Int           *nUentries_ptr,
                                 HYPRE_Int          **Uentries_ptr )
{
   hypre_SStructGraph   *graph   = hypre_SStructMatrixGraph(matrix);
   HYPRE_Int            *split   = hypre_SStructMatrixSplit(matrix, part, var);
   hypre_SStructStencil *stencil = hypre_SStructGraphStencil(graph, part, var);
   HYPRE_Int             entry;
   HYPRE_Int             i;

   HYPRE_Int             nSentries = 0;
   HYPRE_Int            *Sentries  = hypre_SStructMatrixSEntries(matrix);
   HYPRE_Int             nUentries = 0;
   HYPRE_Int            *Uentries  = hypre_SStructMatrixUEntries(matrix);

   for (i = 0; i < nentries; i++)
   {
      entry = entries[i];
      if (entry < hypre_SStructStencilSize(stencil))
      {
         /* stencil entries */
         if (split[entry] > -1)
         {
            Sentries[nSentries] = split[entry];
            nSentries++;
         }
         else
         {
            Uentries[nUentries] = entry;
            nUentries++;
         }
      }
      else
      {
         /* non-stencil entries */
         Uentries[nUentries] = entry;
         nUentries++;
      }
   }

   *nSentries_ptr = nSentries;
   *Sentries_ptr  = Sentries;
   *nUentries_ptr = nUentries;
   *Uentries_ptr  = Uentries;

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * (action > 0): add-to values
 * (action = 0): set values
 * (action < 0): get values
 * (action =-2): get values and zero out
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructMatrixSetValues( HYPRE_SStructMatrix  matrix,
                              HYPRE_Int            part,
                              HYPRE_Int           *index,
                              HYPRE_Int            var,
                              HYPRE_Int            nentries,
                              HYPRE_Int           *entries,
                              HYPRE_Complex       *values,
                              HYPRE_Int            action )
{
   HYPRE_Int             ndim  = hypre_SStructMatrixNDim(matrix);
   hypre_SStructGraph   *graph = hypre_SStructMatrixGraph(matrix);
   hypre_SStructGrid    *grid  = hypre_SStructGraphGrid(graph);
   HYPRE_Int           **nvneighbors = hypre_SStructGridNVNeighbors(grid);
   HYPRE_Int            *Sentries;
   HYPRE_Int            *Uentries;
   HYPRE_Int             nSentries;
   HYPRE_Int             nUentries;
   hypre_SStructPMatrix *pmatrix;
   hypre_Index           cindex;

   hypre_SStructMatrixSplitEntries(matrix, part, var, nentries, entries,
                                   &nSentries, &Sentries,
                                   &nUentries, &Uentries);

   hypre_CopyToCleanIndex(index, ndim, cindex);

   /* S-matrix */
   if (nSentries > 0)
   {
      pmatrix = hypre_SStructMatrixPMatrix(matrix, part);
      hypre_SStructPMatrixSetValues(pmatrix, cindex, var,
                                    nSentries, Sentries, values, action);
      /* put inter-part couplings in UMatrix and zero them out in PMatrix
       * (possibly in ghost zones) */
      if (nvneighbors[part][var] > 0)
      {
         hypre_SStructMatrixSetInterPartValues(matrix, part, cindex, cindex, var,
                                               nSentries, entries, values, action);
      }
   }

   /* U-matrix */
   if (nUentries > 0)
   {
      hypre_SStructUMatrixSetValues(matrix, part, cindex, var,
                                    nUentries, Uentries, values, action);
   }

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * (action > 0): add-to values
 * (action = 0): set values
 * (action < 0): get values
 * (action =-2): get values and zero out
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructMatrixSetBoxValues( HYPRE_SStructMatrix  matrix,
                                 HYPRE_Int            part,
                                 HYPRE_Int           *ilower,
                                 HYPRE_Int           *iupper,
                                 HYPRE_Int            var,
                                 HYPRE_Int            nentries,
                                 HYPRE_Int           *entries,
                                 HYPRE_Complex       *values,
                                 HYPRE_Int            action )
{
   HYPRE_Int                ndim  = hypre_SStructMatrixNDim(matrix);
   hypre_SStructGraph      *graph = hypre_SStructMatrixGraph(matrix);
   hypre_SStructGrid       *grid  = hypre_SStructGraphGrid(graph);
   HYPRE_Int              **nvneighbors = hypre_SStructGridNVNeighbors(grid);
   HYPRE_Int               *Sentries;
   HYPRE_Int               *Uentries;
   HYPRE_Int                nSentries;
   HYPRE_Int                nUentries;
   hypre_SStructPMatrix    *pmatrix;
   hypre_Index              cilower;
   hypre_Index              ciupper;


   hypre_SStructMatrixSplitEntries(matrix, part, var, nentries, entries,
                                   &nSentries, &Sentries,
                                   &nUentries, &Uentries);

   hypre_CopyToCleanIndex(ilower, ndim, cilower);
   hypre_CopyToCleanIndex(iupper, ndim, ciupper);

   /* S-matrix */
   if (nSentries > 0)
   {
      pmatrix = hypre_SStructMatrixPMatrix(matrix, part);
      hypre_SStructPMatrixSetBoxValues(pmatrix, cilower, ciupper, var,
                                       nSentries, Sentries, values, action);

      /* put inter-part couplings in UMatrix and zero them out in PMatrix
       * (possibly in ghost zones) */
      if (nvneighbors[part][var] > 0)
      {
         hypre_SStructMatrixSetInterPartValues(matrix, part, cilower, ciupper, var,
                                               nSentries, entries, values, action);
      }
   }

   /* U-matrix */
   if (nUentries > 0)
   {
      hypre_SStructUMatrixSetBoxValues(matrix, part, cilower, ciupper, var,
                                       nUentries, Uentries, values, action);
   }

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * Put inter-part couplings in UMatrix and zero them out in PMatrix (possibly in
 * ghost zones).  Assumes that all entries are stencil entries.
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_SStructMatrixSetInterPartValues( HYPRE_SStructMatrix  matrix,
                                       HYPRE_Int            part,
                                       hypre_Index          ilower,
                                       hypre_Index          iupper,
                                       HYPRE_Int            var,
                                       HYPRE_Int            nentries,
                                       HYPRE_Int           *entries,
                                       HYPRE_Complex       *values,
                                       HYPRE_Int            action )
{
   HYPRE_Int                ndim  = hypre_SStructMatrixNDim(matrix);
   hypre_SStructGraph      *graph = hypre_SStructMatrixGraph(matrix);
   hypre_SStructGrid       *grid  = hypre_SStructGraphGrid(graph);
   hypre_SStructPMatrix    *pmatrix = hypre_SStructMatrixPMatrix(matrix, part);
   hypre_SStructPGrid      *pgrid = hypre_SStructPMatrixPGrid(pmatrix);
   hypre_BoxArray          *partbnd_boxa = hypre_SStructPGridPBndBoxArray(pgrid, var);
   HYPRE_Int               *smap = hypre_SStructPMatrixSMap(pmatrix, var);
   hypre_SStructVariable    frvartype = hypre_SStructPGridVarType(pgrid, var);
   hypre_SStructStencil    *stencil = hypre_SStructPMatrixStencil(pmatrix, var);
   hypre_Index             *shape = hypre_SStructStencilShape(stencil);
   HYPRE_Int               *vars = hypre_SStructStencilVars(stencil);

   hypre_SStructVariable    tovartype;
   hypre_StructMatrix      *smatrix;
   hypre_Box               *box, *vbox;
   hypre_Box               *ibox0, *ibox1;
   hypre_Box               *tobox, *frbox;

   hypre_Index              stride, loop_size;
   hypre_IndexRef           offset, start;
   hypre_BoxManEntry      **frentries, **toentries;
   hypre_SStructBoxManInfo *frinfo, *toinfo;
   HYPRE_Complex           *tvalues = NULL;
   HYPRE_Int                nfrentries, ntoentries, frpart, topart;
   HYPRE_Int                entry, sentry, ei, fri, toi, vi, mi;

   box   = hypre_BoxCreate(ndim);
   vbox  = hypre_BoxCreate(ndim);
   ibox0 = hypre_BoxCreate(ndim);
   ibox1 = hypre_BoxCreate(ndim);
   tobox = hypre_BoxCreate(ndim);
   frbox = hypre_BoxCreate(ndim);
   hypre_BoxSetExtents(vbox, ilower, iupper);

   hypre_SetIndex(stride, 1);
   for (ei = 0; ei < nentries; ei++)
   {
      entry  = entries[ei];
      sentry = smap[entry];
      offset = shape[entry];
      smatrix = hypre_SStructPMatrixSMatrix(pmatrix, var, vars[entry]);
      tovartype = hypre_SStructPGridVarType(pgrid, vars[entry]);

      /* shift box in the stencil offset direction */
      hypre_BoxSetExtents(box, ilower, iupper);
      hypre_AddIndexes(hypre_BoxIMin(box), offset, ndim, hypre_BoxIMin(box));
      hypre_AddIndexes(hypre_BoxIMax(box), offset, ndim, hypre_BoxIMax(box));

      /* get "to" entries */
      hypre_SStructGridIntersect(grid, part, vars[entry], box, -1,
                                 &toentries, &ntoentries);

      for (toi = 0; toi < ntoentries; toi++)
      {
         hypre_BoxManEntryGetExtents(
            toentries[toi], hypre_BoxIMin(tobox), hypre_BoxIMax(tobox));
         hypre_IntersectBoxes(box, tobox, ibox0);
         if (hypre_BoxVolume(ibox0))
         {
            hypre_SStructBoxManEntryGetPart(toentries[toi], part, &topart);

            /* shift ibox0 back */
            hypre_SubtractIndexes(hypre_BoxIMin(ibox0), offset, ndim,
                                  hypre_BoxIMin(ibox0));
            hypre_SubtractIndexes(hypre_BoxIMax(ibox0), offset, ndim,
                                  hypre_BoxIMax(ibox0));

            /* get "from" entries */
            hypre_SStructGridIntersect(grid, part, var, ibox0, -1,
                                       &frentries, &nfrentries);
            for (fri = 0; fri < nfrentries; fri++)
            {
               /* don't set couplings within the same part unless possibly for
                * cell data (to simplify periodic conditions for users) */
               hypre_SStructBoxManEntryGetPart(frentries[fri], part, &frpart);
               if (topart == frpart)
               {
                  if ( (frvartype != HYPRE_SSTRUCT_VARIABLE_CELL) ||
                       (tovartype != HYPRE_SSTRUCT_VARIABLE_CELL) )
                  {
                     continue;
                  }
                  hypre_BoxManEntryGetInfo(frentries[fri], (void **) &frinfo);
                  hypre_BoxManEntryGetInfo(toentries[toi], (void **) &toinfo);
                  if ( hypre_SStructBoxManInfoType(frinfo) ==
                       hypre_SStructBoxManInfoType(toinfo) )
                  {
                     continue;
                  }
               }

               hypre_BoxManEntryGetExtents(
                  frentries[fri], hypre_BoxIMin(frbox), hypre_BoxIMax(frbox));
               hypre_IntersectBoxes(ibox0, frbox, ibox1);
               if (hypre_BoxVolume(ibox1))
               {
                  tvalues =
                     hypre_TReAlloc(tvalues, HYPRE_Complex, hypre_BoxVolume(ibox1));

                  if (action >= 0)
                  {
                     /* set or add */
                     hypre_AppendBox(ibox1, partbnd_boxa);

                     /* copy values into tvalues */
                     start = hypre_BoxIMin(ibox1);
                     hypre_BoxGetSize(ibox1, loop_size);
                     hypre_BoxLoop2Begin(ndim, loop_size,
                                         ibox1, start, stride, mi,
                                         vbox,  start, stride, vi);
#ifdef HYPRE_USING_OPENMP
#pragma omp parallel for private(HYPRE_BOX_PRIVATE,mi,vi) HYPRE_SMP_SCHEDULE
#endif
                     hypre_BoxLoop2For(mi, vi)
                     {
                        tvalues[mi] = values[ei + vi*nentries];
                     }
                     hypre_BoxLoop2End(mi, vi);

                     /* put values into UMatrix */
                     hypre_SStructUMatrixSetBoxValues(
                        matrix, part, hypre_BoxIMin(ibox1), hypre_BoxIMax(ibox1),
                        var, 1, &entry, tvalues, action);
                     /* zero out values in PMatrix (possibly in ghost) */
                     hypre_StructMatrixClearBoxValues(
                        smatrix, ibox1, 1, &sentry, -1, 1);
                  }
                  else
                  {
                     /* get */

                     /* get values from UMatrix */
                     hypre_SStructUMatrixSetBoxValues(
                        matrix, part, hypre_BoxIMin(ibox1), hypre_BoxIMax(ibox1),
                        var, 1, &entry, tvalues, action);

                     /* copy tvalues into values */
                     start = hypre_BoxIMin(ibox1);
                     hypre_BoxGetSize(ibox1, loop_size);
                     hypre_BoxLoop2Begin(ndim, loop_size,
                                         ibox1, start, stride, mi,
                                         vbox,  start, stride, vi);
#ifdef HYPRE_USING_OPENMP
#pragma omp parallel for private(HYPRE_BOX_PRIVATE,mi,vi) HYPRE_SMP_SCHEDULE
#endif
                     hypre_BoxLoop2For(mi, vi)
                     {
                        values[ei + vi*nentries] = tvalues[mi];
                     }
                     hypre_BoxLoop2End(mi, vi);

                  } /* end if action */
               } /* end if nonzero ibox1 */
            } /* end of "from" boxman entries loop */
            hypre_TFree(frentries);
         } /* end if nonzero ibox0 */
      } /* end of "to" boxman entries loop */
      hypre_TFree(toentries);
   } /* end of entries loop */

   hypre_BoxDestroy(box);
   hypre_BoxDestroy(vbox);
   hypre_BoxDestroy(ibox0);
   hypre_BoxDestroy(ibox1);
   hypre_BoxDestroy(tobox);
   hypre_BoxDestroy(frbox);
   hypre_TFree(tvalues);

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_SStructMatrixToUMatrix
 *
 * Notes (VPM):
 *       1) Use part and var as arguments to this function?
 *       2) We are not converting the whole SStructMatrix, only the
 *          structured part. Change function's name?
 *       3) This converts only A(vi, vi). Need to expand to other variables.
 *--------------------------------------------------------------------------*/

hypre_IJMatrix *
hypre_SStructMatrixToUMatrix( HYPRE_SStructMatrix  matrix )
{
   HYPRE_IJMatrix           ij_A;
   HYPRE_Int                ndim     = hypre_SStructMatrixNDim(matrix);
   HYPRE_Int                nparts   = hypre_SStructMatrixNParts(matrix);
   HYPRE_Int               *Sentries = hypre_SStructMatrixSEntries(matrix);
   MPI_Comm                 comm     = hypre_SStructMatrixComm(matrix);

   hypre_SStructPMatrix    *pmatrix;
   hypre_StructMatrix      *smatrix;
   hypre_StructStencil     *sstencil;
   hypre_StructGrid        *sgrid;
   hypre_BoxArray          *grid_boxes;
   hypre_Box               *box;
   HYPRE_Int               *row_sizes;
   HYPRE_Int               *split;
   HYPRE_Int                nSentries;
   hypre_Index              unit_stride, loop_size;
   hypre_IndexRef           start;
   hypre_IndexRef           ilower, iupper;

   HYPRE_Int                i, m, mi, nrows, ncols, nnzs, max_size;
   HYPRE_Int                istart, iend, jstart, jend;
   HYPRE_Int                part, var, nvars, entry;
   HYPRE_Complex           *values;

   /* Set beggining/end of rows and columns that belong to this process */
   nrows  = hypre_SStructMatrixRanGhlocalSize(matrix);
   ncols  = hypre_SStructMatrixDomGhlocalSize(matrix);
   istart = hypre_SStructMatrixRanGhstartRank(matrix);
   jstart = hypre_SStructMatrixDomGhstartRank(matrix);
   iend   = istart + nrows - 1;
   jend   = jstart + ncols - 1;

   /* Set row sizes */
   max_size = m = 0;
   hypre_SetIndex(unit_stride, 1);
   row_sizes = hypre_CTAlloc(HYPRE_Int, nrows);
   for (part = 0; part < nparts; part++)
   {
      pmatrix = hypre_SStructMatrixPMatrix(matrix, part);
      nvars   = hypre_SStructPMatrixNVars(pmatrix);

      for (var = 0; var < nvars; var++)
      {
         split    = hypre_SStructMatrixSplit(matrix, part, var);
         smatrix  = hypre_SStructPMatrixSMatrix(pmatrix, var, var);
         sstencil = hypre_StructMatrixStencil(smatrix);
         sgrid    = hypre_StructMatrixGrid(smatrix);

         nnzs = 0;
         for (entry = 0; entry < hypre_StructStencilSize(sstencil); entry++)
         {
            if (split[entry] > -1)
            {
               nnzs++;
            }
         }

         grid_boxes = hypre_StructGridBoxes(sgrid);
         hypre_ForBoxI(i, grid_boxes)
         {
            box = hypre_BoxArrayBox(grid_boxes,i);

            start = hypre_BoxIMin(box);
            hypre_BoxGetSize(box, loop_size);
            hypre_BoxLoop1Begin(ndim, loop_size, box, start, unit_stride, mi);
#ifdef HYPRE_USING_OPENMP
#pragma omp parallel for private(HYPRE_BOX_PRIVATE,mi) HYPRE_SMP_SCHEDULE
#endif
            hypre_BoxLoop1For(mi)
            {
               row_sizes[m+mi] = nnzs;
            }
            hypre_BoxLoop1End(mi);

            m += hypre_BoxVolume(box);
         } /* Loop over boxes */

         max_size = hypre_max(max_size, nnzs);
      } /* Loop over vars */
   } /* Loop over parts */

   /* Create and initialize ij_A */
   HYPRE_IJMatrixCreate(comm, istart, iend, jstart, jend, &ij_A);
   HYPRE_IJMatrixSetObjectType(ij_A, HYPRE_PARCSR);
   HYPRE_IJMatrixSetRowSizes(ij_A, (const HYPRE_Int *) row_sizes);
   HYPRE_IJMatrixInitialize(ij_A);

   /* Free/Allocate memory */
   hypre_TFree(row_sizes);
   values = hypre_CTAlloc(HYPRE_Complex, nrows*max_size);

   /* Set entries of ij_A */
   for (part = 0; part < nparts; part++)
   {
      pmatrix = hypre_SStructMatrixPMatrix(matrix, part);
      nvars   = hypre_SStructPMatrixNVars(pmatrix);

      for (var = 0; var < nvars; var++)
      {
         split    = hypre_SStructMatrixSplit(matrix, part, var);
         smatrix  = hypre_SStructPMatrixSMatrix(pmatrix, var, var);
         sstencil = hypre_StructMatrixStencil(smatrix);
         sgrid    = hypre_StructMatrixGrid(smatrix);

         nSentries = 0;
         for (entry = 0; entry < hypre_StructStencilSize(sstencil); entry++)
         {
            if (split[entry] > -1)
            {
               Sentries[nSentries] = split[entry];
               nSentries++;
            }
         }

         grid_boxes = hypre_StructGridBoxes(sgrid);
         hypre_ForBoxI(i, grid_boxes)
         {
            box    = hypre_BoxArrayBox(grid_boxes,i);
            ilower = hypre_BoxIMin(box);
            iupper = hypre_BoxIMax(box);

            /* GET values from this box */
            hypre_SStructPMatrixSetBoxValues(pmatrix, ilower, iupper, var,
                                             nSentries, Sentries, values, -1);

            /* SET values to ij_A */
            hypre_SStructUMatrixSetBoxValuesHelper(matrix, part, ilower, iupper,
                                                   var, nSentries, Sentries,
                                                   values, 0, ij_A);
         } /* Loop over boxes */
      } /* Loop over vars */
   } /* Loop over parts */

   /* Assemble ij_A */
   HYPRE_IJMatrixAssemble(ij_A);

   /* Free memory */
   hypre_TFree(values);

   return (hypre_IJMatrix *) ij_A;
}
