// Command server runs the Delaunay/Voronoi geometry HTTP service.
package main

import (
	"log"
	"os"

	"github.com/example/delaunay/internal/api"
)

func main() {
	port := os.Getenv("PORT")
	if port == "" {
		port = "8080"
	}
	r := api.Router()
	log.Printf("delaunay geometry service listening on :%s", port)
	if err := r.Run(":" + port); err != nil {
		log.Fatalf("server failed: %v", err)
	}
}
