#include <AMReX_Particles.H>
#include <AMReX_AmrCore.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_NeighborParticles.H>

using namespace amrex;

#define AMREX_FI_NSTRUCTREAL BL_SPACEDIM
#define AMREX_FI_NSTRUCTINT 0

namespace {
  using FNborParticleContainer = NeighborParticleContainer<AMREX_FI_NSTRUCTREAL,
                                                       AMREX_FI_NSTRUCTINT,0,0>;
}

extern "C" {

    void amrex_fi_new_nborparticlecontainer (FNborParticleContainer*& particlecontainer,
                                         AmrCore* amrcore)
    {
	particlecontainer = new FNborParticleContainer(amrcore);
    }

    void amrex_fi_delete_nborparticlecontainer (FNborParticleContainer* particlecontainer)
    {
	delete particlecontainer;
    }

    void amrex_fi_get_next_particle_id (int& id)
    {
        id = FNborParticleContainer::ParticleType::NextID();
    }

    void amrex_fi_get_cpu (int& cpu)
    {
        cpu = ParallelDescriptor::MyProc();
    }
    
    void amrex_fi_write_particles(FNborParticleContainer* particlecontainer,
                                  const char* dirname, const char* pname, int is_checkpoint)
    {
        particlecontainer->Checkpoint(dirname, pname, is_checkpoint);
    }

    void amrex_fi_particle_redistributelocal (FNborParticleContainer* particlecontainer,
                                         int lev_min, int lev_max, int ng)
    {
      particlecontainer->RedistributeLocal();
    }

    void amrex_fi_get_particles_mfi(FNborParticleContainer* particlecontainer,
                                    int lev, MFIter* mfi, Real*& dp, Long& np)
    {
        const int grid = mfi->index();
        const int tile = mfi->LocalTileIndex();
        auto& particle_level = particlecontainer->GetParticles(lev);
        auto search = particle_level.find(std::make_pair(grid, tile));
        if (search != particle_level.end()) {
            auto& particle_tile = search->second;
            np = particle_tile.numParticles();
            if (np > 0) {
                auto& aos = particle_tile.GetArrayOfStructs();
                dp = aos.data();
            } else {
                dp = nullptr;
            }
        } else {
            np = 0;
            dp = nullptr;
        }
    }
    
    void amrex_fi_add_particle_mfi(FNborParticleContainer* particlecontainer,
                                   int lev, MFIter* mfi, FNborParticleContainer::ParticleType* p)
    {
        const int grid = mfi->index();
        const int tile = mfi->LocalTileIndex();
        auto& particle_level = particlecontainer->GetParticles(lev);
        auto& particle_tile  = particle_level[std::make_pair(grid, tile)];
        particle_tile.push_back(*p);
    }
    
    void amrex_fi_num_particles_mfi(FNborParticleContainer* particlecontainer,
                                    int lev, MFIter* mfi, Long& np)
    {
        const int grid = mfi->index();
        const int tile = mfi->LocalTileIndex();
        auto& particle_level = particlecontainer->GetParticles(lev);
        auto search = particle_level.find(std::make_pair(grid, tile));
        if (search != particle_level.end()) {
            auto& particle_tile = search->second;
            np = particle_tile.numParticles();
        } else {
            np = 0;
        }
    }

    void amrex_fi_get_particles_i(FNborParticleContainer* particlecontainer,
                                  int lev, int grid, int tile, Real*& dp, Long& np)
    {
        auto& particle_level = particlecontainer->GetParticles(lev);
        auto search = particle_level.find(std::make_pair(grid, tile));
        if (search != particle_level.end()) {
            auto& particle_tile = search->second;
            np = particle_tile.numParticles();
            if (np > 0) {
                auto& aos = particle_tile.GetArrayOfStructs();
                dp = aos.data();
            } else {
                dp = nullptr;
            }
        } else {
            np = 0;
            dp = nullptr;
        }
    }
    
    void amrex_fi_add_particle_i(FNborParticleContainer* particlecontainer,
                                 int lev, int grid, int tile, FNborParticleContainer::ParticleType* p)
    {
        auto& particle_level = particlecontainer->GetParticles(lev);
        auto& particle_tile  = particle_level[std::make_pair(grid, tile)];
        particle_tile.push_back(*p);
    }
    
    void amrex_fi_num_particles_i(FNborParticleContainer* particlecontainer,
                                  int lev, int grid, int tile, Long& np)
    {
        auto& particle_level = particlecontainer->GetParticles(lev);
        auto search = particle_level.find(std::make_pair(grid, tile));
        if (search != particle_level.end()) {            
            auto& particle_tile = search->second;
            np = particle_tile.numParticles();
        } else {
            np = 0;
        }
    }

  void amrex_fi_nborparticle_regrid(FNborParticleContainer* particlecontainer, DistributionMapping& dm, BoxArray& ba, int lev)
  {
    particlecontainer->Regrid(dm, ba, lev);
  }

  void amrex_fi_nborparticle_fillNeighbors (FNborParticleContainer* Particlecontainer)
  {
    particlecontainer->fillNeighbors();
  }

  void amrex_fi_nborparticle_sumNeighbors(FNborParticleContainer* particlecontainer, int real_start_comp, int real_num_comp, int int_start_comp, int int_num_comp)
  {
    particlecontainer->sumNeighbors(real_start_comp, real_num_comp, int_start_comp, int_num_comp);
  }

  void amrex_fi_nborparticle_updateNeighbors(FNborParticleContainer* particlecontainer)
  {
    particlecontainer->updateNeighbors();
  }

  void amrex_fi_nborparticle_clearNeighbors(FNborParticleContainer* particlecontainer)
  {
    particlecontainer->clearNeighbors();
  }

  void amrex_fi_nborparticle_RedistributeLocal(FNborParticleContainer* particlecontainer)
  {
    particlecontainer->RedistributeLocal();
  }

  void amrex_fi_nborparticle_buildNeighborList(FNborParticleContainer* particlecontainer, CheckPair&& check_pair)
  {
    bool sort=false;
    particlecontainer->buildNeighborList(check_pair, sort);
  }

  void amrex_fi_nborparticle_printNeighborList(FNborParticleContainer* particlecontainer)
  {
    particlecontainer->printNeighborList();
  }

  void amrex_fi_nborparticle_GetNeighbors(FNborParticleContainer* particlecontainer, ParticleVector& neighbors, int lev, int grid, int tile)
  {
    return particlecontainer->GetNeighbors(lev, grid, tile);
  }

}
