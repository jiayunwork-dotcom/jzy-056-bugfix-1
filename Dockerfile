# syntax=docker/dockerfile:1

# ---- build stage: pinned to Go 1.22 ----
FROM golang:1.22-bookworm AS build

WORKDIR /src

# Resolve dependencies first for better layer caching.
COPY go.mod go.sum ./
RUN go mod download

# Build the static service binary.
COPY . .
ENV CGO_ENABLED=0 GOOS=linux
RUN go build -trimpath -ldflags="-s -w" -o /out/delaunay ./cmd/server

# ---- runtime stage: minimal, no shell needed ----
FROM gcr.io/distroless/static-debian12:nonroot

WORKDIR /
COPY --from=build /out/delaunay /delaunay

ENV PORT=8080
EXPOSE 8080
USER nonroot:nonroot

ENTRYPOINT ["/delaunay"]
